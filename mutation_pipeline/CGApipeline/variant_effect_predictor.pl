#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

http://www.ensembl.org/info/docs/tools/vep/script/index.html

Version 83

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use CGI qw/:standard/;
use FindBin qw($RealBin);
use lib $RealBin;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VEP qw(
    parse_line
    vf_to_consequences
    validate_vf
    convert_to_vcf
    get_all_consequences
    get_slice
    build_full_cache
    read_cache_info
    get_version_data
    get_time
    debug
    @OUTPUT_COLS
    @VCF_COLS
    @EXTRA_HEADERS
    %COL_DESCS
    @REG_FEAT_TYPES
    %FILTER_SHORTCUTS
    @PICK_ORDER
);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

# global vars
my $VERSION = '83';

my %ts_tv = (
  'A/G' => 'Ts',
  'G/A' => 'Ts',
  'C/T' => 'Ts',
  'T/C' => 'Ts',
  'A/C' => 'Tv',
  'C/A' => 'Tv',
  'G/T' => 'Tv',
  'T/G' => 'Tv',
  'C/G' => 'Tv',
  'G/C' => 'Tv',
  'A/T' => 'Tv',
  'T/A' => 'Tv',
);

my %colour_keys = (
  'polyphen' => {
    'unknown' => 'blue',
    'benign' => 'green',
    'possibly damaging' => 'orange',
    'probably damaging' => 'red',
  },
  'sift' => {
    'tolerated' => 'green',
    'deleterious' => 'red',
  },
  
  # copied from COLOUR.ini in web code via browser to check colours
  'consequences' => {
    'intergenic_variant'                => 'gray',
    'intron_variant'                    => '#02599c',
    'upstream_gene_variant'             => '#a2b5cd',
    'downstream_gene_variant'           => '#a2b5cd',
    '5_prime_utr_variant'               => '#7ac5cd',
    '3_prime_utr_variant'               => '#7ac5cd',
    'splice_region_variant'             => '#ff7f50',
    'splice_donor_variant'              => '#ff7f50',
    'splice_acceptor_variant'           => '#ff7f50',
    'frameshift_variant'                => '#ff69b4',
    'transcript_ablation'               => '#ff0000',
    'transcript_amplification'          => '#ff69b4',
    'inframe_insertion'                 => '#ff69b4',
    'inframe_deletion'                  => '#ff69b4',
    'synonymous_variant'                => '#76ee00',
    'stop_retained_variant'             => '#76ee00',
    'missense_variant'                  => '#ffd700',
    'initiator_codon_variant'           => '#ffd700',
    'stop_gained'                       => '#ff0000',
    'stop_lost'                         => '#ff0000',
    'mature_mirna_variant'              => '#458b00',
    'non_coding_exon_variant'           => '#32cd32',
    'nc_transcript_variant'             => '#32cd32',
    'incomplete_terminal_codon_variant' => '#ff00ff',
    'nmd_transcript_variant'            => '#ff4500',
    'coding_sequence_variant'           => '#458b00',
    'tfbs_ablation'                     => 'brown',
    'tfbs_amplification'                => 'brown',
    'tf_binding_site_variant'           => 'brown',
    'regulatory_region_variant'         => 'brown',
    'regulatory_region_ablation'        => 'brown',
    'regulatory_region_amplification'   => 'brown',
  },
);

# set output autoflush for progress bars
$| = 1;

# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# run with leaktrace
# use Test::LeakTrace;
# leaktrace { &main($config); undef $config; } -verbose;

# this is the main sub-routine - it needs the configured $config hash
sub main {
  my $config = shift;
    
  debug("Starting...") unless defined $config->{quiet};
    
  # this is for counting seconds
  $config->{start_time} = time();
  $config->{last_time} = time();
    
  # this is for stats
  $config->{stats}->{start_time} = get_time();
    
  my $tr_cache = {};
  my $rf_cache = {};
    
  # create a hash to hold slices so we don't get the same one twice
  my %slice_cache = ();
    
  my @vfs;    
  my ($vf_count, $total_vf_count);
  my $in_file_handle = $config->{in_file_handle};
    
  # initialize line number in config
  $config->{line_number} = 0;
    
  # read the file
  while(<$in_file_handle>) {
    
    # split again to avoid Windows character nonsense
    foreach my $line(split /\r|(?>\v|\x0D\x0A)/) {
      
      chomp($line);
      
      next unless $line =~ /\w+/;
        
      $config->{line_number}++;
        
      # header line?
      if($line =~ /^\#/) {
            
        # retain header lines if we are outputting VCF
        if(defined($config->{vcf})) {
          push @{$config->{headers}}, $line;
        }
            
        # line with sample labels in VCF
        if(defined($config->{individual}) && /^#CHROM/) {
          my @split = split(/\s+/, $line);
                
          # no individuals
          die("ERROR: No individual data found in VCF\n") if scalar @split <= 9;
                
          # get individual column indices
          my %ind_cols = map {$split[$_] => $_} (9..$#split);
                
          # all?
          if(scalar @{$config->{individual}} == 1 && $config->{individual}->[0] =~ /^all$/i) {
            $config->{ind_cols} = \%ind_cols;
          }
          else {
            my %new_ind_cols;
                    
            # check we have specified individual(s)
            foreach my $ind(@{$config->{individual}}) {
              die("ERROR: Individual named \"$ind\" not found in VCF\n") unless defined $ind_cols{$ind};
              $new_ind_cols{$ind} = $ind_cols{$ind};
            }
                    
            $config->{ind_cols} = \%new_ind_cols;
          }
        }
            
        next;
      }
        
      # strip off nasty characters
      $line =~ s/\s+$//g;
        
      # configure output file
      $config->{out_file_handle} ||= &get_out_file_handle($config);
        
      # some lines (pileup) may actually parse out into more than one variant
      foreach my $vf(@{&parse_line($config, $line)}) {
            
        $vf->{_line} = $line;
        $vf->{_line_number} = $config->{line_number};
            
        # now get the slice
        if(!defined($vf->{slice})) {
          my $slice;
                
          # don't get slices if we're using cache
          # we can steal them from transcript objects later
          if((!defined($config->{cache}) && !defined($config->{whole_genome})) || defined($config->{check_ref}) || defined($config->{convert})) {
                    
            # check if we have fetched this slice already
            if(defined $slice_cache{$vf->{chr}}) {
              $slice = $slice_cache{$vf->{chr}};
            }
                    
            # if not create a new one
            else {
                        
              $slice = &get_slice($config, $vf->{chr});
                        
              # if failed, warn and skip this line
              if(!defined($slice)) {
                warn("WARNING: Could not fetch slice named ".$vf->{chr}." on line ".$config->{line_number}."\n") unless defined $config->{quiet};
                next;
              }    
                        
              # store the hash
              $slice_cache{$vf->{chr}} = $slice;
            }
          }
                
          $vf->{slice} = $slice;
        }
            
        # validate the VF
        next unless validate_vf($config, $vf) || defined($config->{dont_skip});
            
        # make a name if one doesn't exist
        $vf->{variation_name} ||= ($vf->{original_chr} || $vf->{chr}).'_'.$vf->{start}.'_'.($vf->{allele_string} || $vf->{class_SO_term});
            
        # jump out to convert here
        if(defined($config->{convert})) {
          &convert_vf($config, $vf);
          next;
        }
            
        if(defined $config->{whole_genome}) {
          push @vfs, $vf;
          $vf_count++;
          $total_vf_count++;
                
          if($vf_count == $config->{buffer_size}) {
            debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
            
            $config->{stats}->{out_count} += print_line($config, $_) foreach @{get_all_consequences($config, \@vfs)};
                    
            # calculate stats
            my $total_rate = sprintf("%.0f vars/sec", $total_vf_count / ((time() - $config->{start_time}) || 1));
            my $rate = sprintf("%.0f vars/sec", $vf_count / ((time() - $config->{last_time}) || 1));
            $config->{last_time} = time();
                    
            debug("Processed $total_vf_count total variants ($rate, $total_rate total)") unless defined($config->{quiet});
                    
            @vfs = ();
            $vf_count = 0;
          }
        }
        else {
          $config->{stats}->{out_count} += print_line($config, $_) foreach @{vf_to_consequences($config, $vf)};
          $vf_count++;
          $total_vf_count++;
          debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
        }
      }
    }
  }
    
  # if in whole-genome mode, finish off the rest of the buffer
  if(defined $config->{whole_genome} && scalar @vfs) {
    debug("Read $vf_count variants into buffer") unless defined($config->{quiet});
        
    $config->{stats}->{out_count} += print_line($config, $_) foreach @{get_all_consequences($config, \@vfs)};
        
    # calculate stats
    my $total_rate = sprintf("%.0f vars/sec", $total_vf_count / ((time() - $config->{start_time}) || 1));
    my $rate = sprintf("%.0f vars/sec", $vf_count / ((time() - $config->{last_time}) || 1));
    $config->{last_time} = time();
        
    debug("Processed $total_vf_count total variants ($rate, $total_rate total)") unless defined($config->{quiet});
  }
    
  debug($config->{stats}->{filter_count}, "/$total_vf_count variants remain after filtering") if (defined($config->{check_frequency})) && !defined($config->{quiet});
    
  debug("Executed ", defined($Bio::EnsEMBL::DBSQL::StatementHandle::count_queries) ? $Bio::EnsEMBL::DBSQL::StatementHandle::count_queries : 'unknown number of', " SQL statements") if defined($config->{count_queries}) && !defined($config->{quiet});
    
  # finalise run-time stats
  $config->{stats}->{var_count} = $total_vf_count;
  $config->{stats}->{end_time} = get_time();
  $config->{stats}->{run_time} = time() - $config->{start_time};
    
  # write stats
  unless(defined($config->{no_stats})) {
    summarise_stats($config);
    debug("Wrote stats summary to ".$config->{stats_file}) unless defined($config->{quiet});
  }
    
  # tell user about any warnings
  if($config->{warning_count}) {
    debug("See ".$config->{warning_file}." for details of ".$config->{warning_count}." warnings") unless defined($config->{quiet});
  }
    
  # close HTML output
  if(defined($config->{html}) && defined($config->{html_file_handle})) {
    my $fh = $config->{html_file_handle};
    print $fh "</tbody><tfoot><tr>".$config->{_th}."</tr></tfoot></table><p>&nbsp;</p></div></html></body>\n</html>\n";
    $fh->close;
  }
    
  if(defined($config->{solr})) {
    my $fh = $config->{out_file_handle};
    print $fh "</add>\n";
  }
    
  # tabix?
  if(defined($config->{tabix})) {
    debug("Compressing and indexing output") unless defined($config->{quiet});
      
    # check sorting
    open SORT, $config->{output_file};
    my $prev_pos = 0;
    my $prev_chr = 0;
    my $is_sorted = 0;
    my %seen_chrs;
      
    while(<SORT>) {
      if(!/^\#/) {
        $is_sorted = 1;
        my @data = split /\s+/, $_;
        if(
          ($data[0] eq $prev_chr && $data[1] < $prev_pos) ||
          ($data[0] ne $prev_chr && defined($seen_chrs{$data[0]}))
        ) {
          $is_sorted = 0;
          last;
        }
          
        $prev_pos = $data[1];
        $prev_chr = $data[0];
        $seen_chrs{$data[0]} = 1;
      }
    }
      
    close SORT;
      
    unless($is_sorted) {
      system('grep "^#" '.$config->{output_file}.' > '.$config->{output_file}.'.sorted');
      system('grep -v "^#" '.$config->{output_file}.' | sort -k1,1 -k2,2n >> '.$config->{output_file}.'.sorted');
      system('mv '.$config->{output_file}.'.sorted '.$config->{output_file});
    }
      
    system("bgzip -f ".$config->{output_file}) == 0 or warn "WARNING: failed to generated bgzipped file for tabix\n$?\n";
    rename($config->{output_file}.".gz", $config->{output_file});
    system("tabix -p vcf -f ".$config->{output_file}) == 0 or warn "WARNING: failed to generated tabix index\n$?\n";
  }
    
  debug("Finished!") unless defined $config->{quiet};
}

# sets up configuration hash that is used throughout the script
sub configure {
    my $args = shift;
    
    my $config = {};
    
    my @ARGV_copy = @ARGV;
    $config->{stats}->{options} = \@ARGV_copy;
    
    GetOptions(
        $config,
        'help',                    # displays help message
        
        # input options,
        'config=s',                # config file name
        'input_file|i=s',          # input file name
        'format=s',                # input file format
        
        # DB options
        'species=s',               # species e.g. human, homo_sapiens
        'registry=s',              # registry file
        'host=s',                  # database host
        'port=s',                  # database port
        'user|u=s',                  # database user name
        'password=s',              # database password
        'db_version=i',            # Ensembl database version to use e.g. 62
        'assembly|a=s',            # assembly version to use
        'genomes',                 # automatically sets DB params for e!Genomes
        'refseq',                  # use otherfeatures RefSeq DB instead of Ensembl
        'merged',                  # use merged cache
        'all_refseq',              # report consequences on all transcripts in RefSeq cache, includes CCDS, EST etc
        'gencode_basic',           # limit to using just GenCode basic transcript set
       
        'is_multispecies=i',       # '1' for a multispecies database (e.g protists_euglenozoa1_collection_core_29_82_1)
        # runtime options
        'minimal',                 # convert input alleles to minimal representation
        'most_severe',             # only return most severe consequence
        'summary',                 # only return one line per variation with all consquence types
        'per_gene',                # only return most severe per gene
        'pick',                    # used defined criteria to return most severe line
        'flag_pick',               # like --pick but just adds a flag to picked line
        'pick_allele',             # like --pick but chooses one con per allele
        'flag_pick_allele',        # like --flag_pick but flags one con per allele
        'pick_order=s',            # define the order of categories used by the --*pick* flags
        'buffer_size=i',           # number of variations to read in before analysis
        'chunk_size=s',            # size in bases of "chunks" used in internal hash structure
        'failed=i',                # include failed variations when finding existing
        'no_whole_genome',         # disables now default whole-genome mode
        'whole_genome',            # proxy for whole genome mode - now just warns user
        'gp',                      # read coords from GP part of INFO column in VCF (probably only relevant to 1KG)
        'chr=s',                   # analyse only these chromosomes, e.g. 1-5,10,MT
        'check_ref',               # check supplied reference allele against DB
        'check_existing',          # find existing co-located variations
        'check_svs',               # find overlapping structural variations
        'check_alleles',           # only attribute co-located if alleles are the same
        'check_frequency',         # enable frequency checking
        'gmaf',                    # add global MAF of existing var
        'maf_1kg',                 # add 1KG MAFs of existing vars
        'maf_esp',                 # add ESP MAFs of existing vars
        'maf_exac',                # add ExAC MAFs of existing vars
        'old_maf',                 # report 1KG/ESP MAFs in the old way (no allele, always < 0.5)
        'pubmed',                  # add Pubmed IDs for publications that cite existing vars
        'freq_filter=s',           # exclude or include
        'freq_freq=f',             # frequency to filter on
        'freq_gt_lt=s',            # gt or lt (greater than or less than)
        'freq_pop=s',              # population to filter on
        'filter_common',           # shortcut to MAF filtering
        'allow_non_variant',       # allow non-variant VCF lines through
        'process_ref_homs',        # force processing of individuals with homozygous ref genotype
        'individual=s',            # give results by genotype for individuals
        'phased',                  # force VCF genotypes to be interpreted as phased
        'fork=i',                  # fork into N processes
        'dont_skip',               # don't skip vars that fail validation
        
        # verbosity options
        'verbose|v',               # print out a bit more info while running
        'quiet',                   # print nothing to STDOUT (unless using -o stdout)
        'no_progress',             # don't display progress bars
        
        # output options
        'everything|e',            # switch on EVERYTHING :-)
        'output_file|o=s',         # output file name
        'tabix',                   # bgzip and tabix-index output
        'html',                    # generate an HTML version of output
        'stats_file|sf=s',         # stats file name
        'stats_text',              # write stats as text
        'warning_file=s',          # file to write warnings to
        'no_stats',                # don't write stats file
        'force_overwrite',         # force overwrite of output file if already exists
        'terms|t=s',               # consequence terms to use e.g. NCBI, SO
        'coding_only',             # only return results for consequences in coding regions
        'canonical',               # indicates if transcript is canonical
        'tsl',                     # output transcript support level
        'appris',                  # output APPRIS transcript annotation
        'ccds',                    # output CCDS identifer
        'xref_refseq',             # output refseq mrna xref
        'uniprot',                 # output Uniprot identifiers (includes UniParc)
        'protein',                 # add e! protein ID to extra column
        'biotype',                 # add biotype of transcript to output
        'hgnc',                    # add HGNC gene ID to extra column
        'symbol',                  # add gene symbol (e.g. HGNC)
        'gene_phenotype',          # indicate if genes are phenotype-associated
        'hgvs',                    # add HGVS names to extra column
        'shift_hgvs=i',            # disable/enable 3-prime shifting of HGVS indels to comply with standard
        'sift=s',                  # SIFT predictions
        'polyphen=s',              # PolyPhen predictions
        'humdiv',                  # use humDiv instead of humVar for PolyPhen
        'condel=s',                # Condel predictions
        'variant_class',           # get SO variant type
        'regulatory',              # enable regulatory stuff
        'cell_type=s' => ($config->{cell_type} ||= []),             # filter cell types for regfeats
        'convert=s',               # convert input to another format (doesn't run VEP)
        'no_intergenic',           # don't print out INTERGENIC consequences
        'gvf',                     # produce gvf output
        'vcf',                     # produce vcf output
        'solr',                    # produce XML output for Solr
        'json',                    # produce JSON document output
        'vcf_info_field=s',        # allow user to change VCF info field name
        'keep_csq',                # don't nuke existing CSQ fields in VCF
        'keep_ann',                # synonym for keep_csq
        'no_consequences',         # don't calculate consequences
        'lrg',                     # enable LRG-based features
        'fields=s',                # define your own output fields
        'domains',                 # output overlapping protein features
        'numbers',                 # include exon and intron numbers
        'total_length',            # give total length alongside positions e.g. 14/203
        'allele_number',           # indicate allele by number to avoid confusion with VCF conversions
        'no_escape',               # don't percent-escape HGVS strings
        
        # cache stuff
        'database',                # must specify this to use DB now
        'cache',                   # use cache
        'cache_version=i',         # specify a different cache version
        'write_cache',             # enables writing to the cache
        'show_cache_info',         # print cache info and quit
        'build=s',                 # builds cache from DB from scratch; arg is either all (all top-level seqs) or a list of chrs
        'build_test',              # disable some slow start-up stuff for speed when testing
        'build_parts=s',           # choose which bits of the cache to build (t=transcript, v=variants, r=regfeats)
        'build_range=s',           # for testing, give a coord range. Probably best when using one chrom only e.g. --build 21
        'no_adaptor_cache',        # don't write adaptor cache
        'strip',                   # strips adaptors etc from objects before caching them
        'rebuild=s',               # rebuilds cache by reading in existing then redumping - probably don't need to use this any more
        'dir=s',                   # dir where cache is found (defaults to $HOME/.vep/)
        'dir_cache=s',             # specific directory for cache
        'dir_plugins=s',           # specific directory for plugins
        'cache_region_size=i',     # size of region in bases for each cache file
        'no_slice_cache',          # tell API not to cache features on slice
        'standalone',              # standalone renamed offline
        'offline',                 # offline mode uses minimal set of modules installed in same dir, no DB connection
        'skip_db_check',           # don't compare DB parameters with cached
        'compress=s',              # by default we use zcat to decompress; user may want to specify gzcat or "gzip -dc"
        'custom=s' => ($config->{custom} ||= []), # specify custom tabixed bgzipped file with annotation
        'tmpdir=s',                # tmp dir used for BigWig retrieval
        'plugin=s' => ($config->{plugin} ||= []), # specify a method in a module in the plugins directory
        'safe',                    # die if plugins don't compile or spit warnings
        'fasta=s',                 # file or dir containing FASTA files with reference sequence
        'freq_file=s',             # file containing freqs to add to cache build
        'freq_vcf=s' => ($config->{freq_vcf} ||= []), # VCF file containing freqs
        'sereal',                  # user Sereal instead of Storable for the cache
        
        # debug
        'cluck',                   # these two need some mods to Bio::EnsEMBL::DBSQL::StatementHandle to work. Clucks callback trace and SQL
        'count_queries',           # counts SQL queries executed
        'admin',                   # allows me to build off public hosts
        'debug',                   # print out debug info
    ) or die "ERROR: Failed to parse command-line flags\n";
    
    # print usage message if requested or no args supplied
    if(defined($config->{help}) || !$args) {
        &usage;
        exit(0);
    }
    
    # config file?
    if(defined $config->{config}) {
        read_config_from_file($config, $config->{config});
    }
    
    # dir is where the cache and plugins live
    my $default_dir = join '/', ($ENV{'HOME'}, '.vep');
    $config->{dir_plugins} ||= ($config->{dir} ? $config->{dir}.'/Plugins' : $default_dir.'/Plugins');
    $config->{dir} ||= $config->{dir_cache} || $default_dir;

    # ini file?
    my $ini_file = $config->{dir}.'/vep.ini';
    
    if(-e $ini_file) {
        read_config_from_file($config, $ini_file);
    }
    
    # everything?
    if(defined($config->{everything})) {
        my %everything = (
            sift           => 'b',
            polyphen       => 'b',
            ccds           => 1,
            hgvs           => 1,
            symbol         => 1,
            numbers        => 1,
            domains        => 1,
            regulatory     => 1,
            canonical      => 1,
            protein        => 1,
            biotype        => 1,
            gmaf           => 1,
            maf_1kg        => 1,
            maf_esp        => 1,
            maf_exac       => 1,
            pubmed         => 1,
            uniprot        => 1,
            tsl            => 1,
            appris         => 1,
            variant_class  => 1,
            gene_phenotype => 1,
        );
        
        $config->{$_} = $everything{$_} for keys %everything;
    }
    
    # subroutine to check for illegal flags or combinations
    check_flags($config);
    
    # connection settings for Ensembl Genomes
    if($config->{genomes}) {
        $config->{host} ||= 'mysql-eg-publicsql.ebi.ac.uk';
        $config->{port} ||= 4157;
    }
    
    # connection settings for main Ensembl
    else {
        $config->{species} ||= "homo_sapiens";
        $config->{host}    ||= 'ensembldb.ensembl.org';
        $config->{port}    ||= 3306;
    }
    
    # refseq or core?
    if(defined($config->{refseq})) {
        $config->{core_type} = 'otherfeatures';
    }
    else {
        $config->{core_type} = 'core';
    }
    
    # turn on rest for json
    if(defined($config->{json})) {
      $config->{rest} = 1;
      
      eval q{ use JSON; };
      
      if($@) {
        die("ERROR: Could not load required JSON module\n");
      }
    }
    
    if(defined($config->{vcf})) {
      $config->{keep_csq} = 1 if defined($config->{keep_ann});
      $config->{$_} = 1 for qw(symbol biotype numbers);
    }
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    # individual(s) specified?
    if(defined($config->{individual})) {
        $config->{individual} = [split /\,/, $config->{individual}];
        
        # force allow_non_variant
        $config->{allow_non_variant} = 1;
    }
    
    # regulatory has to be on for cell_type
    if(defined($config->{cell_type}) && scalar(@{$config->{cell_type}})) {
        $config->{regulatory} = 1;
        $config->{cell_type} = [map {split /\,/, $_} @{$config->{cell_type}}];
    }
    else {
      delete $config->{cell_type};
    }
    
    # summarise options if verbose
    if(defined $config->{verbose}) {
        my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
        print $header;
        
        my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
        
        foreach my $key(sort keys %$config) {
            next if ref($config->{$key}) eq 'ARRAY' && scalar @{$config->{$key}} == 0;
            print $key.(' ' x (($max_length - length($key)) + 4)).(ref($config->{$key}) eq 'ARRAY' ? join "\t", @{$config->{$key}} : $config->{$key})."\n";
        }
        
        print "\n".("-" x 20)."\n\n";
    }
    
    # set defaults
    $config->{user}              ||= 'anonymous';
    $config->{buffer_size}       ||= 5000;
    $config->{chunk_size}        ||= '50kb';
    $config->{output_file}       ||= "variant_effect_output.txt";
    $config->{stats_file}        ||= $config->{output_file}."_summary.".(defined($config->{stats_text}) ? 'txt' : 'html');
    $config->{tmpdir}            ||= '/tmp';
    $config->{format}            ||= 'guess';
    $config->{terms}             ||= 'SO';
    $config->{cache_region_size} ||= 1000000;
    $config->{compress}          ||= 'gzip -dc';
    $config->{polyphen_analysis}   = defined($config->{humdiv}) ? 'humdiv' : 'humvar';
    $config->{vcf_info_field}    ||= 'CSQ';
    
    # shift HGVS?
    if(defined($config->{shift_hgvs})) {
      use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
      use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
      no warnings 'once';
      $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $config->{shift_hgvs};
      no warnings 'once';
      $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = $config->{shift_hgvs};
    }
    
    # frequency filtering
    if(defined($config->{filter_common})) {
        $config->{check_frequency} = 1;
        
        # set defaults
        $config->{freq_freq}   ||= 0.01;
        $config->{freq_filter} ||= 'exclude';
        $config->{freq_pop}    ||= '1KG_ALL';
        $config->{freq_gt_lt}  ||= 'gt';
    }
    
    if(defined($config->{check_frequency})) {
        foreach my $flag(qw(freq_freq freq_filter freq_pop freq_gt_lt)) {
            die "ERROR: To use --check_frequency you must also specify flag --$flag\n" unless defined $config->{$flag};
        }
        
        # need to set check_existing
        $config->{check_existing} = 1;
    }
    
    foreach my $flag(qw(check_existing check_alleles gmaf maf_1kg maf_esp maf_exac pubmed)) {
      $config->{check_existing} = 1 if defined $config->{$flag};
    }
    
    # warn users still using whole_genome flag
    if(defined($config->{whole_genome})) {
        debug("INFO: Whole-genome mode is now the default run-mode for the script. To disable it, use --no_whole_genome") unless defined($config->{quiet});
    }
    
    $config->{whole_genome}      = 1 unless defined $config->{no_whole_genome};
    $config->{failed}            = 0 unless defined $config->{failed};
    $config->{chunk_size}        =~ s/mb?/000000/i;
    $config->{chunk_size}        =~ s/kb?/000/i;
    $config->{cache_region_size} =~ s/mb?/000000/i;
    $config->{cache_region_size} =~ s/kb?/000/i;
    
    # cluck and display executed SQL?
    if (defined($config->{cluck})) {
        no warnings 'once';
        $Bio::EnsEMBL::DBSQL::StatementHandle::cluck = 1;
    }
    
    # write_cache needs cache
    $config->{cache} = 1 if defined $config->{write_cache};
    $config->{cache} = 1 if defined $config->{offline};
    $config->{cache} = 1 if defined $config->{show_cache_info};
    
    # no_slice_cache and whole_genome have to be on to use cache
    if(defined($config->{cache})) {
        $config->{no_slice_cache} = 1;
        $config->{whole_genome} = 1;
        $config->{strip} = 1;
    }
    
    $config->{build} = $config->{rebuild} if defined($config->{rebuild});
    
    # force options for full build
    if(defined($config->{build})) {
        $config->{symbol} = 1;
        $config->{no_slice_cache} = 1;
        $config->{cache} = 1;
        $config->{strip} = 1;
        $config->{write_cache} = 1;
        $config->{cell_type} = [1] if defined($config->{regulatory});
        $config->{build_parts} ||= 'tvr';
        
        die("ERROR: --build_parts [tvr] must have at least one of t (transcripts), v (variants), r (regfeats)\n") unless $config->{build_parts} =~ /^[tvr]+$/;
    } 
    
    # connect to databases
    $config->{reg} = &connect_to_dbs($config);
    
    # setup cache dir etc
    setup_cache($config) if defined($config->{cache});
    
    if(defined($config->{show_cache_info})) {
      my $v = get_version_data($config);
      print "$_\t".$v->{$_}."\n" for keys %$v;
      exit(0);
    }
    
    # setup FASTA file
    if(defined($config->{fasta})) {
      # spoof a coordinate system
      $config->{coord_system} = Bio::EnsEMBL::CoordSystem->new(
        -NAME => 'chromosome',
        -RANK => 1,
      );

      $config->{fasta_db} = setup_fasta(
        -FASTA => $config->{fasta},
        -ASSEMBLY => $config->{assembly},
        -OFFLINE => $config->{offline},
      ) if defined($config->{fasta});
    }
    
    # setup custom files
    setup_custom($config) if defined($config->{custom});
    
    # setup forking
    setup_forking($config) if defined($config->{fork});
    
    # offline needs cache, can't use HGVS
    if(defined($config->{offline})) {
        
        die("ERROR: Cannot generate HGVS coordinates in offline mode without a FASTA file (see --fasta)\n") if defined($config->{hgvs}) && !defined($config->{fasta});
        die("ERROR: Cannot use HGVS as input in offline mode\n") if $config->{format} eq 'hgvs';
        die("ERROR: Cannot use variant identifiers as input in offline mode\n") if $config->{format} eq 'id';
        die("ERROR: Cannot do frequency filtering in offline mode\n") if defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
        die("ERROR: Cannot retrieve overlapping structural variants in offline mode\n") if defined($config->{check_sv});
        die("ERROR: Cannot check reference sequences without a FASTA file (see --fasta)\n") if defined($config->{check_ref}) && !defined($config->{fasta});
        die("ERROR: Cannot map to LRGs in offline mode\n") if defined($config->{lrg});
    }
    
    # suppress warnings that the FeatureAdpators spit if using no_slice_cache
    Bio::EnsEMBL::Utils::Exception::verbose(1999) if defined($config->{no_slice_cache});
   
    # we configure plugins here because they can sometimes switch on the 
    # regulatory config option
    configure_plugins($config);
    
    # include regulatory modules if requested
    if(defined($config->{regulatory})) {
        # do the use statements here so that users don't have to have the
        # funcgen API installed to use the rest of the script
        eval q{
            use Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor;
            use Bio::EnsEMBL::Funcgen::MotifFeature;
            use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
            use Bio::EnsEMBL::Funcgen::BindingMatrix;
        };
        
        if($@) {
            die("ERROR: Ensembl Funcgen API must be installed to use --regulatory or plugins that deal with regulatory features\n$@");
        }
    }
    
    # user defined custom output fields
    if(defined($config->{fields})) {
        $config->{fields} = [split ',', $config->{fields}];
        debug("Output fields redefined (".scalar @{$config->{fields}}." defined)") unless defined($config->{quiet});
        $config->{fields_redefined} = 1;
    }
    $config->{fields} ||= \@OUTPUT_COLS;
    
    # get adaptors (don't get them in offline mode)
    unless(defined($config->{offline})) {
        
        &get_adaptors($config);
        
        # reg adaptors (only fetches if not retrieved from cache already)
        &get_reg_adaptors($config) if defined($config->{regulatory});
    }
    
    # check regulatory available
    if(defined($config->{regulatory}) && defined($config->{cache}) && !defined($config->{write_cache}) && !(defined($config->{cache_regulatory}) || defined($config->{cache_cell_types}))) {
        
        # --everything this option is implicit so don't die
        if(defined($config->{everything})) {
            delete $config->{regulatory};
        }
        else {
            die("ERROR: --regulatory is not available for this species");
        }
    }
    
    # check cell types
    if(defined($config->{cell_type}) && scalar @{$config->{cell_type}} && !defined($config->{build})) {
        my $cls = '';
        
        if(defined($config->{cache})) {
            $cls = $config->{cache_cell_types};
        }
        else {
            my $cta = $config->{RegulatoryFeature_adaptor}->db->get_CellTypeAdaptor();
            $cls = join ",", map {$_->name} @{$cta->fetch_all};
        }
        
        foreach my $cl(@{$config->{cell_type}}) {
            die "ERROR: cell type $cl not recognised; available cell types are:\n$cls\n" unless $cls =~ /(^|,)$cl(,|$)/;
        }
    }
    
    # check SIFT/PolyPhen available?
    foreach my $tool(grep {defined($config->{$_})} qw(sift polyphen)) {
      my $vd = get_version_data($config);
        
      unless(defined($vd->{$tool}) || defined($config->{'cache_'.$tool.'_version'})) {
            
        # --everything this option is implicit so don't die
        if(defined($config->{everything})) {
          delete $config->{$tool};
        }
        else {
          die("ERROR: --$tool is not available for this species or cache");
        }
      }
    }
    
    # get terminal width for progress bars
    unless(defined($config->{quiet}) || defined($config->{no_progress})) {
        my $width;
        
        # module may not be installed
        eval q{
            use Term::ReadKey;
        };
        
        if(!$@) {
            my ($w, $h);
            
            # module may be installed, but e.g.
            eval {
                ($w, $h) = GetTerminalSize();
            };
            
            $width = $w if defined $w;
        }
        
        $width ||= 60;
        $width -= 12;
        $config->{terminal_width} = $width;
    }
    
    # jump out to build cache if requested
    if(defined($config->{build})) {
        
        if($config->{host} =~ /^(ensembl|useast)db\.ensembl\.org$/ && !defined($config->{admin})) {
            die("ERROR: Cannot build cache using public database server ", $config->{host}, "\n");
        }
        
        # get VCF freqs
        if(defined($config->{freq_vcf})) {
          my @new;
          
          foreach my $vcf_conf(@{$config->{freq_vcf}}) {
            my ($freq_file, @opts_and_pops) = split /\,/, $vcf_conf;

            my @opts = grep {/\=/} @opts_and_pops;
            my @file_pops = grep {!/\=/} @opts_and_pops;

            my %opts = (
              file => $freq_file,
              pops => \@file_pops,
            );

            foreach my $opt(@opts) {
              my ($k, $v) = split('=', $opt);
              $opts{$k} = $v;
            }

            # create prefixed pop names if given
            my $prefix = $opts{prefix} || '';
            $prefix .= '_' if $prefix && $prefix !~ /\_$/;
            @{$opts{prefixed_pops}} = map {s/\_$//; $_} map {$prefix.$_} @{$opts{pops}};
            
            push @new, \%opts;
            
            push @{$config->{'freq_file_pops'}}, @{$opts{prefixed_pops}};
          }
          
          $config->{freq_vcf} = \@new;
        }
        
        if(defined($config->{'freq_file'})) {
            my ($freq_file, @file_pops) = split /\,/, $config->{'freq_file'};
            debug("Loading extra frequencies from $freq_file") unless defined($config->{quiet});
            
            open IN, $freq_file or die "ERROR: Could not open frequencies file $freq_file\n";
            while(<IN>) {
                chomp;
                
                if(m/^(\w+)(\s)/) {
                  my $id = $1;
                  s/^$1$2//;
                  tr/\t/ /;
                  $config->{'freqs'}->{$id} = $_;
                }
            }
            close IN;
            
            # add pops to $config
            push @{$config->{'freq_file_pops'}}, @file_pops;

            push @{$config->{'just_file_pops'}}, @file_pops;
        }
        
        # build the cache
        debug("Building cache for ".$config->{species}) unless defined($config->{quiet});
        build_full_cache($config);
        
        # exit script
        debug("Finished building cache") unless defined($config->{quiet});
        exit(0);
    }
    
    # warn user DB will be used for SIFT/PolyPhen/HGVS/frequency/LRG
    if(defined($config->{cache})) {
        
        # these two def depend on DB
        foreach my $param(grep {defined $config->{$_}} qw(hgvs lrg check_sv check_ref)) {
            debug("INFO: Database will be accessed when using --$param") unless defined($config->{quiet}) or ($param =~ /hgvs|check_ref/ and defined($config->{fasta}));
        }
        
        debug("INFO: Database will be accessed when using --check_frequency with population ".$config->{freq_pop}) if !defined($config->{quiet}) and defined($config->{check_frequency}) && $config->{freq_pop} !~ /1kg.*(all|afr|amr|asn|eur)/i;
        
        # as does using HGVS or IDs as input
        debug("INFO: Database will be accessed when using --format ", $config->{format}) if ($config->{format} eq 'id' || $config->{format} eq 'hgvs') && !defined($config->{quiet});
    }
    
    # get list of chrs if supplied
    if(defined($config->{chr})) {
        my %chrs;
        
        foreach my $val(split /\,/, $config->{chr}) {
            my @nnn = split /\-/, $val;
            
            foreach my $chr($nnn[0]..$nnn[-1]) {
                $chrs{$chr} = 1;
            }
        }
        
        $config->{chr} = \%chrs;
    }
    
    # get input file handle
    $config->{in_file_handle} = &get_in_file_handle($config);
    
    return $config;
}

# reads config from a file
sub read_config_from_file {
    my $config = shift;
    my $file = shift;
    
    open CONFIG, $file or die "ERROR: Could not open config file \"$file\"\n";
    
    while(<CONFIG>) {
        next if /^\#/;
        
        # preserve spaces between quotes
        s/([\"\'].*)(\s)(.*[\"\'])/$1\_\_\_SPACE\_\_\_$3/g;
        
        my @split = split /\s+|\=/;
        my $key = shift @split;
        $key =~ s/^\-//g;
        
        # restore spaces
        s/\_\_\_SPACE\_\_\_/ /g for @split;
        
        # remove quotes
        s/[\"\']//g for @split;
        
        if(defined($config->{$key}) && ref($config->{$key}) eq 'ARRAY') {
            push @{$config->{$key}}, @split;
        }
        else {
            $config->{$key} ||= $split[0];
        }
    }
    
    close CONFIG;
    
    # force quiet if outputting to STDOUT
    if(defined($config->{output_file}) && $config->{output_file} =~ /stdout/i) {
        delete $config->{verbose} if defined($config->{verbose});
        $config->{quiet} = 1;
    }
    
    debug("Read configuration from $file") unless defined($config->{quiet});
}

# configures custom VEP plugins
sub configure_plugins {

    my $config = shift;
    
    $config->{plugins} = [];
    
    if (my @plugins = @{ $config->{plugin} }) {

        # add the Plugins directory onto @INC

        unshift @INC, $config->{dir_plugins};

        for my $plugin (@plugins) {

            # parse out the module name and parameters

            my ($module, @params) = split /,/, $plugin;

            # check we can use the module
            
            eval qq{
                use $module;
            };
            if ($@) {
                my $msg = "Failed to compile plugin $module: $@\n";
                die($msg) if defined($config->{safe});
                debug($msg) unless defined($config->{quiet});
                next;
            }
            
            # now check we can instantiate it, passing any parameters to the constructor
            
            my $instance;
            
            eval {
                $instance = $module->new($config, @params);
            };
            if ($@) {
                my $msg = "Failed to instantiate plugin $module: $@\n";
                die($msg) if defined($config->{safe});
                debug($msg) unless defined($config->{quiet});
                next;
            }

            # check that the versions match
            
            #my $plugin_version;
            #
            #if ($instance->can('version')) {
            #    $plugin_version = $instance->version;
            #}
            #
            #my $version_ok = 1;
            #
            #if ($plugin_version) {
            #    my ($plugin_major, $plugin_minor, $plugin_maintenance) = split /\./, $plugin_version;
            #    my ($major, $minor, $maintenance) = split /\./, $VERSION;
            #
            #    if ($plugin_major != $major) {
            #        debug("Warning: plugin $plugin version ($plugin_version) does not match the current VEP version ($VERSION)") unless defined($config->{quiet});
            #        $version_ok = 0;
            #    }
            #}
            #else {
            #    debug("Warning: plugin $plugin does not define a version number") unless defined($config->{quiet});
            #    $version_ok = 0;
            #}
            #
            #debug("You may experience unexpected behaviour with this plugin") unless defined($config->{quiet}) || $version_ok;

            # check that it implements all necessary methods
            
            for my $required(qw(run get_header_info check_feature_type check_variant_feature_type)) {
                unless ($instance->can($required)) {
                    my $msg = "Plugin $module doesn't implement a required method '$required', does it inherit from BaseVepPlugin?\n";
                    die($msg) if defined($config->{safe});
                    debug($msg) unless defined($config->{quiet});
                    next;
                }
            }
           
            # all's good, so save the instance in our list of plugins
            
            push @{ $config->{plugins} }, $instance;
            
            debug("Loaded plugin: $module") unless defined($config->{quiet}); 

            # for convenience, check if the plugin wants regulatory stuff and turn on the config option if so
            
            if (grep { $_ =~ /motif|regulatory/i } @{ $instance->feature_types }) {
                debug("Fetching regulatory features for plugin: $module") unless defined($config->{quiet});
                $config->{regulatory} = 1;
            }
        }
    }
} 

# connects to DBs (not done in offline mode)
sub connect_to_dbs {
    my $config = shift;
    
    # get registry
    my $reg = 'Bio::EnsEMBL::Registry';
    
    unless(defined($config->{offline})) {
        # load DB options from registry file if given
        if(defined($config->{registry})) {
            debug("Loading DB config from registry file ", $config->{registry}) unless defined($config->{quiet});
            $reg->load_all(
                $config->{registry},
                $config->{verbose},
                undef,
                $config->{no_slice_cache}
            );
        }
        
        # otherwise manually connect to DB server
        else {
	  if($config->{is_multispecies}==1){
             $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
             );
           }else {
       	     $reg->load_registry_from_db(
                -host       => $config->{host},
                -user       => $config->{user},
                -pass       => $config->{password},
                -port       => $config->{port},
                -db_version => $config->{db_version},
                -species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
                -verbose    => $config->{verbose},
                -no_cache   => $config->{no_slice_cache},
            );
	  }
        }
        
        eval { $reg->set_reconnect_when_lost() };
        
        # get meta container adaptors to check version
        my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
        my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
        
        if(defined($config->{verbose})) {
            
            if($core_mca && $var_mca) {
                debug(
                    "Connected to core version ", $core_mca->get_schema_version, " database ",
                    "and variation version ", $var_mca->get_schema_version, " database"
                );
            }
        }
        
        # get assembly version
        if($core_mca) {
          my ($highest_cs) = @{$core_mca->db->get_CoordSystemAdaptor->fetch_all()};
          my $assembly = $highest_cs->version();
          
          die(
            "ERROR: Assembly version specified by --assembly (".$config->{assembly}.
            ") and assembly version in coord_system table (".$assembly.") do not match\n".
            (
              $config->{host} eq 'ensembldb.ensembl.org' ?
              "\nIf using human GRCh37 add \"--port 3337\"".
              " to use the GRCh37 database, or --offline to avoid database connection entirely\n" :
              ''
            )
          ) if defined($config->{assembly}) && $config->{assembly} ne $assembly;
          
          $config->{assembly} = $assembly;
          
          if(!defined($config->{assembly})) {
            die("ERROR: No assembly version specified, use --assembly [version] or check the coord_system table in your core database\n");
          }
        }
    }
    
    return $reg;
}

# get adaptors from DB
sub get_adaptors {
    my $config = shift;
    
    die "ERROR: No registry" unless defined $config->{reg};
    
    # try fetching a slice adaptor
    eval {
        $config->{sa}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'slice');
    };
    
    if($@) {
        if($@ =~ /not find internal name for species/) {
            my %register = %Bio::EnsEMBL::Registry::registry_register;
            my @species_list = sort keys %{$register{_SPECIES}};
            
            debug("ERROR: Could not find database for species ".$config->{species});
            
            if(scalar @species_list) {
                debug("List of valid species for this server:\n\n".(join "\n", @species_list));
            }
            
            die("\nExiting\n");
        }
        
        else {
            die $@;
        }
    }
    
    # get the remaining core adaptors
    $config->{ga}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'gene');
    $config->{ta}  = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'transcript');
    $config->{mca} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'metacontainer');
    $config->{csa} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'coordsystem');
    $config->{tra} = $config->{reg}->get_adaptor($config->{species}, $config->{core_type}, 'translation');
    
    # get variation adaptors
    $config->{vfa}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variationfeature');
    $config->{svfa}  = $config->{reg}->get_adaptor($config->{species}, 'variation', 'structuralvariationfeature');
    $config->{tva}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'transcriptvariation');
    $config->{pfpma} = $config->{reg}->get_adaptor($config->{species}, 'variation', 'proteinfunctionpredictionmatrix');
    $config->{va}    = $config->{reg}->get_adaptor($config->{species}, 'variation', 'variation');
    $config->{pfa}   = $config->{reg}->get_adaptor($config->{species}, 'variation', 'phenotypefeature');
    
    # get fake ones for species with no var DB
    if(!defined($config->{vfa})) {
        $config->{vfa}  = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($config->{species});
        $config->{svfa} = Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor->new_fake($config->{species});
        $config->{tva}  = Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor->new_fake($config->{species});
    }
    
    # cache schema version
    $config->{mca}->get_schema_version if defined $config->{mca};
    
    # check we got slice adaptor - can't continue without a core DB
    die("ERROR: Could not connect to core database\n") unless defined $config->{sa};
}

sub check_flags() {
  my $config = shift;
  
  my @invalid = (
    ['quiet', 'verbose'],
    ['refseq', 'gencode_basic'],
    ['refseq', 'merged'],
    ['merged', 'database'],
    ['database', 'offline'],
    ['database', 'cache'],
  );

  foreach my $combo(@invalid) {
    die "ERROR: Can't use these flags together: ".join(", ", map {'--'.$_} @$combo)."\n" if scalar(grep {defined($config->{$_})} @$combo) == scalar @$combo;
  }
  
  # required
  # die "ERROR: --all_refseq requires using either --refseq or --merged\n" if defined($config->{all_refseq}) && !defined($config->{refseq}) && !defined($config->{merged});
  
  # check for deprecated flags
  die "ERROR: --hgnc has been replaced by --symbol\n" if defined($config->{hgnc});
  die "ERROR: --standalone replaced by --offline\n" if(defined $config->{standalone});
  
  # check one of database/cache/offline/build
  if(!grep {defined($config->{$_})} qw(database cache offline build convert show_cache_info)) {
    die qq{
IMPORTANT INFORMATION:

The VEP can read gene data from either a local cache or local/remote databases.

Using a cache is the fastest and most efficient way to use the VEP. The
included INSTALL.pl script can be used to fetch and set up cache files from the
Ensembl FTP server. Simply run "perl INSTALL.pl" and follow the instructions, or
see the documentation pages listed below.

If you have already set up a cache, use "--cache" or "--offline" to use it.

It is possible to use the public databases hosted at ensembldb.ensembl.org, but
this is slower than using the cache and concurrent and/or long running VEP jobs
can put strain on the Ensembl servers, limiting availability to other users.

To enable using databases, add the flag "--database".

Documentation
Installer: http://www.ensembl.org/info/docs/tools/vep/vep_script.html#installer
Cache: http://www.ensembl.org/info/docs/tools/vep/script/index.html#cache

    }
  };
  
  foreach my $flag(qw(maf_1kg maf_esp maf_exac pubmed)) {
    die("ERROR: \-\-$flag can only be used with --cache or --offline")
      if defined($config->{$flag}) && !(defined($config->{cache}) || defined($config->{offline}));
  }
  
  # output term
  if(defined $config->{terms}) {
    die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
    if($config->{terms} =~ /ensembl|display/i) {
      $config->{terms} = 'display';
    }
    else {
      $config->{terms} = uc($config->{terms});
    }
  }
  
  # check file format
  if(defined $config->{format}) {
    die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /^(pileup|vcf|guess|hgvs|ensembl|id)$/i;
  }
  
  # check convert format
  if(defined $config->{convert}) {
    die "ERROR: Unrecognised output format for conversion specified \"".$config->{convert}."\"\n" unless $config->{convert} =~ /vcf|ensembl|pileup|hgvs/i;
    
    # disable stats
    $config->{no_stats} = 1;
  }
  
  # check nsSNP tools
  foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
    die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
    
    die "ERROR: $tool functionality is now available as a VEP Plugin - see http://www.ensembl.org/info/docs/variation/vep/vep_script.html#plugins\n" if $tool eq 'Condel';
  }
  
  # output format has to be VCF for tabix
  die "ERROR: Output must be vcf (--vcf) to use --tabix\n" if defined($config->{tabix}) && !defined($config->{vcf});
  
  # can't use more than one of most_severe, pick, per_gene, summary
  my $total_sev_opts = 0;
  map {$total_sev_opts++ if defined($config->{$_})} qw(most_severe pick pick_allele flag_pick per_gene summary);
  die "ERROR: Can't use more than one of --most_severe, --pick, --per_gene, --summary\n" if $total_sev_opts > 1;
  
  # can't use a whole bunch of options with most_severe
  if(defined($config->{most_severe})) {
    foreach my $flag(qw(no_intergenic protein symbol sift polyphen coding_only ccds canonical xref_refseq numbers domains summary pick_order)) {
      die "ERROR: --most_severe is not compatible with --$flag\n" if defined($config->{$flag});
    }
  }
  
  # can't use a whole bunch of options with summary
  if(defined($config->{summary})) {
    foreach my $flag(qw(no_intergenic protein symbol sift polyphen coding_only ccds canonical xref_refseq numbers domains most_severe pick_order)) {
      die "ERROR: --summary is not compatible with --$flag\n" if defined($config->{$flag});
    }
  }
  
  # check pick order
  if(defined($config->{pick_order})) {
    $config->{pick_order} = [split ',', $config->{pick_order}];
    
    my %valid = map {$_ => 1} @PICK_ORDER;
    my $valid_str = join(", ", @PICK_ORDER);
    
    foreach my $cat(@{$config->{pick_order}}) {
      die("ERROR: $cat is not a valid category for --pick_order, valid categories are:\n$valid_str\n") unless defined($valid{$cat});
    }
  }
}

sub setup_cache() {
  my $config = shift;   
  
  # complete dir with species name and db_version
  my $species_dir_name = defined($config->{offline}) ? $config->{species} : ($config->{reg}->get_alias($config->{species}) || $config->{species});
  $species_dir_name .= '_refseq' if defined($config->{refseq});
  $species_dir_name .= '_merged' if defined($config->{merged});
  
  # add species dir name
  $config->{dir} .= '/'.$species_dir_name;
  
  # check whats in here to match to assembly if given
  die("ERROR: Cache directory ", $config->{dir}, " not found\n") if !-e $config->{dir} && !defined($config->{write_cache});
  
  my $cache_version = $config->{cache_version} || $config->{db_version} || $config->{reg}->software_version;
  
  opendir DIR, $config->{dir};
  my @dir_contents = grep {!/^\./} readdir DIR;
  closedir DIR;
  
  # writing to cache?
  if(defined($config->{write_cache})) {
    
    if(!defined($config->{assembly})) {
      die("ERROR: No assembly specified, or assembly version not found in core DB meta table\n");
    }
    
    else {
      $config->{dir} .= '/'.$cache_version.'_'.$config->{assembly};
      
      if(-e $config->{dir}) {
        debug("INFO: Existing cache directory ", $config->{dir}, " found - contents may be overwritten") unless defined($config->{quiet});
      }
      
      else {
        debug("INFO: Cache directory ", $config->{dir}, " not found - it will be created") unless defined($config->{quiet});
      }
    }
  }
  
  # just reading from cache
  else {
    my @matched_contents = grep {/^$cache_version/} @dir_contents;
    
    # no matched entries, cache not installed
    if(scalar @matched_contents == 0) {
      die("ERROR: No cache found for $species_dir_name, version $cache_version\n");
    }
    
    # only 1 entry, can assume this is OK
    elsif(scalar @matched_contents == 1) {
      $config->{dir} .= '/'.$matched_contents[0];
      
      # is there an assembly here?
      if($matched_contents[0] =~ /\d+\_(.+)/) {
        my $matched_assembly = $1;
        
        if(defined($config->{assembly}) && $config->{assembly} ne $matched_assembly) {
          die(
            "ERROR: Cache assembly version ($matched_assembly) ".
            "and database or selected assembly version (".$config->{assembly}.
            ") do not match\n".
            (
              $config->{host} eq 'ensembldb.ensembl.org' ?
              "\nIf using human GRCh37 add \"--port 3337\"".
              " to use the GRCh37 database, or --offline to avoid database connection entirely\n" :
              ''
            )
          );
        }
        
        $config->{assembly} ||= $matched_assembly;
      }
    }
    
    # did user specify assembly version?
    elsif(!defined($config->{assembly})) {
      my $possibles = join(", ", map {s/^$cache_version\_//; $_} @matched_contents);
      die("ERROR: Multiple assemblies found for cache version $cache_version ($possibles) - specify one using --assembly [assembly]\n");
    }
    
    # add cache version and assembly
    else {
      $config->{dir} .= '/'.$cache_version.'_'.$config->{assembly};
      die("ERROR: No cache found for ".$config->{species}.", version $cache_version, assembly ".$config->{assembly}."\n") unless -e $config->{dir};
    }
  }
  
  # read cache info
  if(read_cache_info($config)) {
    debug("Read existing cache info") unless defined $config->{quiet};
  }
  
  # check assembly matches
  if(defined($config->{assembly}) && defined($config->{cache_assembly}) && $config->{assembly} ne $config->{cache_assembly}) {
    die("ERROR: Mismatch in assembly versions from config (".$config->{assembly}.") and cache info.txt file (".$config->{cache_assembly}.")\n");
  }
  
  # check if there's a FASTA file in there
  if(!defined($config->{fasta})) {
    opendir CACHE, $config->{dir};
    my ($fa) = grep {/\.fa(\.gz)?$/} readdir CACHE;
    
    if(defined $fa) {
      $config->{fasta} = $config->{dir}.'/'.$fa;
      debug("Auto-detected FASTA file in cache directory") unless defined $config->{quiet};
    }
  }
  
  # disable HGVS if no FASTA file found and it was switched on by --everything
  if(
    defined($config->{hgvs}) &&
    defined($config->{offline}) &&
    !defined($config->{fasta}) &&
    defined($config->{everything})
  ) {
    debug("INFO: Disabling --hgvs; using --offline and no FASTA file found\n");
    delete $config->{hgvs};
  }
  
  # check if any disabled options are in use
  # these are set in the cache info file
  if(defined($config->{cache_disabled})) {
    my @arr = ref($config->{cache_disabled} eq 'ARRAY') ? @{$config->{cache_disabled}} : ($config->{cache_disabled});
    
    if(my ($disabled) = grep {defined($config->{$_})} @arr) {
      die("ERROR: Unable to use --".$disabled." with this cache\n");
    }
  }

  # enable sereal
  if(defined($config->{cache_serialiser_type}) && $config->{cache_serialiser_type} eq 'sereal') {
    $config->{sereal} = 1;
  }

  if(defined($config->{sereal})) {
    eval q{ use Sereal; };
    die("ERROR: Could not use Sereal perl module; perhaps you forgot to install it?\n$@") if $@;
  }
}


sub setup_custom {
  my $config = shift;
  
  # check custom annotations
  for my $i(0..$#{$config->{custom}}) {
    my $custom = $config->{custom}->[$i];
    
    my ($filepath, $shortname, $format, $type, $coords, @fields) = split /\,/, $custom;
    $type ||= 'exact';
    $format ||= 'bed';
    $coords ||= 0;
    
    # check type
    die "ERROR: Type $type for custom annotation file $filepath is not allowed (must be one of \"exact\", \"overlap\")\n" unless $type =~ /exact|overlap/;
    
    # check format
    die "ERROR: Format $format for custom annotation file $filepath is not allowed (must be one of \"bed\", \"vcf\", \"gtf\", \"gff\", \"bigwig\")\n" unless $format =~ /bed|vcf|gff|gtf|bigwig/;
    
    # bigwig format
    if($format eq 'bigwig') {
      # check for bigWigToWig
      die "ERROR: bigWigToWig does not seem to be in your path - this is required to use bigwig format custom annotations\n" unless `which bigWigToWig 2>&1` =~ /bigWigToWig$/;
    }
    
    else {
      # check for tabix
      die "ERROR: tabix does not seem to be in your path - this is required to use custom annotations\n" unless `which tabix 2>&1` =~ /tabix$/;
    
      # remote files?
      my $filename = (split /\//, $filepath)[-1];
    
      if($filepath =~ /tp\:\/\//) {
        if(!-e $filename.'.tbi') {
          my $remote_test = `tabix -f $filepath 1:1-1 2>&1`;
          
          if($remote_test =~ /get_local_version/) {
            debug("Downloaded tabix index file for remote annotation file $filepath") unless defined($config->{quiet});
          }
          else {
            die "$remote_test\nERROR: Could not find file or index file for remote annotation file $filepath\n";
          }
        }
      }
  
      # check files exist
      else {
        die "ERROR: Custom annotation file $filepath not found\n" unless -e $filepath;
        die "ERROR: Tabix index file $filepath\.tbi not found - perhaps you need to create it first?\n" unless -e $filepath.'.tbi';
      }
    }
    
    $config->{custom}->[$i] = {
      'file'   => $filepath,
      'name'   => $shortname || 'CUSTOM'.($i + 1),
      'type'   => $type,
      'format' => $format,
      'coords' => $coords,
      'fields' => \@fields
    };
  }
}

sub setup_forking {
  my $config = shift;
  
  if($config->{fork} == 0) {
    delete $config->{fork};
    return;
  }

  die "ERROR: Fork number must be greater than 1\n" if $config->{fork} <= 1;
  
  # check we can use MIME::Base64
  eval q{ use MIME::Base64; };
  
  if($@) {
    debug("WARNING: Unable to load MIME::Base64, forking disabled") unless defined($config->{quiet});
    delete $config->{fork};
  }
  else {
    
    # try a practice fork
    my $pid = fork;
    
    if(!defined($pid)) {
      debug("WARNING: Fork test failed, forking disabled") unless defined($config->{quiet});
      delete $config->{fork};
    }
    elsif($pid) {
      waitpid($pid, 0);
    }
    elsif($pid == 0) {
      exit(0);
    }
  }
}

# gets regulatory adaptors
sub get_reg_adaptors {
    my $config = shift;

    foreach my $type(@REG_FEAT_TYPES) {
        next if defined($config->{$type.'_adaptor'});
        
        my $adaptor = $config->{reg}->get_adaptor($config->{species}, 'funcgen', $type);
        if(defined($adaptor)) {
            $config->{$type.'_adaptor'} = $adaptor;
        }
        else {
            die("ERROR: --regulatory is not available for this species");
        }
    }
}

# gets file handle for input
sub get_in_file_handle {
    my $config = shift;

    # define the filehandle to read input from
    my $in_file_handle = new FileHandle;
    
    if(defined($config->{input_file})) {
        
        # check defined input file exists
        die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
        
        if(-B $config->{input_file}){
            $in_file_handle->open($config->{compress}." ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
        else {
            $in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
        }
    }
    
    # no file specified - try to read data off command line
    else {
        $in_file_handle = 'STDIN';
        debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
    }
    
    return $in_file_handle;
}

# gets file handle for output and adds header
sub get_out_file_handle {
    my $config = shift;
    
    # define filehandle to write to
    my $out_file_handle = new FileHandle;
    
    # check if file exists
    if(-e $config->{output_file} && !defined($config->{force_overwrite})) {
        die("ERROR: Output file ", $config->{output_file}, " already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
    }
    
    if($config->{output_file} =~ /stdout/i) {
        $out_file_handle = *STDOUT;
    }
    #elsif(defined($config->{tabix})) {
    #    $out_file_handle->open(" | bgzip -c > ".$config->{output_file});
    #}
    else {
        $out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
    }
    
    # get stats file handle
    unless(defined $config->{no_stats}) {
      
      # do same for stats file
      if(-e $config->{stats_file} && !defined($config->{force_overwrite})) {
          die("ERROR: Stats file ", $config->{stats_file}, " already exists. Specify a different output file with --stats_file or overwrite existing file with --force_overwrite\n");
      }
      
      die("ERROR: Stats file name ", $config->{stats_file}, " doesn't end in \".htm\" or \".html\" - some browsers may not be able to open this file\n") unless $config->{stats_file} =~ /htm(l)?$/ || defined($config->{stats_text});
      my $stats_file_handle = new FileHandle;
      $stats_file_handle->open(">".$config->{stats_file}) or die("ERROR: Could not write to stats file ", $config->{stats_file}, "\n");
      $config->{stats_file_handle} = $stats_file_handle;
    }
    
    # HTML output?
    my $html_file_handle;
    
    if(defined($config->{html})) {
      if(-e $config->{output_file}.'.html' && !defined($config->{force_overwrite})) {
          die("ERROR: HTML file ", $config->{output_file}, ".html already exists. Specify a different output file with --output_file or overwrite existing file with --force_overwrite\n");
      }
      
      $html_file_handle = new FileHandle;
      $html_file_handle->open(">".$config->{output_file}.'.html') or die("ERROR: Could not write to HTML file ", $config->{output_file}, ".html\n");
      $config->{html_file_handle} = $html_file_handle;
      
      print $html_file_handle html_head();
    }
    
    # define headers for a VCF file
    my @vcf_headers = (
        '#CHROM',
        'POS',
        'ID',
        'REF',
        'ALT',
        'QUAL',
        'FILTER',
        'INFO'
    );
    
    # file conversion, don't want to add normal headers
    if(defined($config->{convert})) {
        # header for VCF
        if($config->{convert} =~ /vcf/i) {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
        }
        
        return $out_file_handle;
    }
    
    # GVF output, no header
    elsif(defined($config->{gvf})) {
        return $out_file_handle;
    }
    
    # VCF format
    elsif(defined($config->{vcf})) {
        
        # create an info string for the VCF header        
        my (@new_headers, @vcf_info_strings);
        
        my $vcf_version_string = sprintf(
          "##VEP=v%i cache=%s db=%s",
          $VERSION,
          $config->{cache} ? $config->{dir} : '.',
          $config->{offline} ? '.' : (
            $config->{mca} ? $config->{mca}->dbc->dbname.'@'.$config->{mca}->dbc->host : '.'
          )
        );
        
        my $version_data = get_version_data($config);
        $vcf_version_string .= ' '.$_.'='.$version_data->{$_} for keys %$version_data;
        
        push @vcf_info_strings, $vcf_version_string;
        
        # if the user has defined the fields themselves, we don't need to worry
        if(defined $config->{fields_redefined}) {
            @new_headers = @{$config->{fields}};
        }
        else {
            
            # now we have to reconfigure headers to comply with the snpEff ANN standard
            my %vcf_cols = map {$_ => 1} @VCF_COLS;
            
            @new_headers = @VCF_COLS;
            
            push @new_headers, (
                
                grep {!$vcf_cols{$_}}
                
                # get default headers, minus variation name and location (already encoded in VCF)
                grep {
                    $_ ne 'Uploaded_variation' and
                    $_ ne 'Location' and
                    $_ ne 'Extra'
                } @{$config->{fields}},
                
                # get extra headers
                map {@{$_->{cols}}}
                grep {defined $config->{$_->{flag}}}
                @EXTRA_HEADERS
            );
            
            # plugin headers
            foreach my $plugin_header(@{get_plugin_headers($config)}) {
                my ($key, $value) = @$plugin_header;
                push @vcf_info_strings, sprintf('##%s=%s', $key, $value);
                push @new_headers, $key;
            }
            
            # redefine the main headers list in config
            $config->{fields} = \@new_headers;
        }
        
        # add the newly defined headers as a header to the VCF
        my $string = join '|', @{$config->{fields}};
        push @vcf_info_strings, '##INFO=<ID='.$config->{vcf_info_field}.',Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '.$string.'">';
        
        # add custom headers
        foreach my $custom(@{$config->{custom}}) {
            push @vcf_info_strings, '##INFO=<ID='.$custom->{name}.',Number=.,Type=String,Description="'.$custom->{file}.' ('.$custom->{type}.')">';
            
            foreach my $field(@{$custom->{fields}}) {
              push @vcf_info_strings, '##INFO=<ID='.$custom->{name}.'_'.$field.',Number=.,Type=String,Description="'.$field.' field from '.$custom->{file}.'">';
            }
        }
        
        # if this is already a VCF file, we need to add our new headers in the right place
        if(defined($config->{headers})) {
            
            # nuke existing CSQ header unless we are keeping it
            unless(defined($config->{keep_csq})) {
              my $vcf_field = $config->{vcf_info_field};
              @{$config->{headers}} = grep {!/^##INFO=<ID=$vcf_field,/} @{$config->{headers}};
            }
            
            for my $i(0..$#{$config->{headers}}) {
                if($config->{headers}->[$i] =~ /^\#CHROM\s+POS\s+ID/) {
                    splice(@{$config->{headers}}, $i, 0, @vcf_info_strings);
                }
            }
            
            print $out_file_handle join "\n", @{$config->{headers}};
            print $out_file_handle "\n";
            
            if(defined($config->{html})) {
                my @tmp = @{$config->{headers}};
                my @cols = split /\s+/, pop @tmp;
                print $html_file_handle join "\n", @tmp;
                print $html_file_handle "\n";
                print $html_file_handle html_table_headers($config, \@cols);
            }
        }
        
        else {
            print $out_file_handle "##fileformat=VCFv4.0\n";
            print $out_file_handle join "\n", @vcf_info_strings;
            print $out_file_handle "\n";
            print $out_file_handle join "\t", @vcf_headers;
            print $out_file_handle "\n";
            
            if(defined($config->{html})) {
                print $html_file_handle "##fileformat=VCFv4.0\n";
                print $html_file_handle join "\n", @vcf_info_strings;
                print $html_file_handle "\n";
                print $html_file_handle html_table_headers($config, \@vcf_headers);
            }
        }
        
        return $out_file_handle;
    }
    
    elsif(defined($config->{solr}) || defined($config->{json})) {
        my @new_headers;
        
        # if the user has defined the fields themselves, we don't need to worry
        if(!defined $config->{fields_redefined}) {
            @new_headers = (
                
                # get default headers
                grep {$_ ne 'Extra'} @{$config->{fields}},
                
                # get extra headers
                map {@{$_->{cols}}}
                grep {defined $config->{$_->{flag}}}
                @EXTRA_HEADERS
            );
            
            # plugin headers
            push @new_headers, map {$_->[0]} @{get_plugin_headers($config)};
            
            # redefine the main headers list in config
            $config->{fields} = \@new_headers;
        }
        
        print $out_file_handle "<add>\n" if defined($config->{solr});
        return $out_file_handle;
    }
    
    # make header
    my $time = &get_time;
    my $db_string = "## Connected to ".$config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host if defined $config->{mca};
    $db_string .= ($db_string ? "\n" : "")."## Using cache in ".$config->{dir} if defined($config->{cache});
    my $version_string =
        "Using API version ".$config->{reg}->software_version.
        ", DB version ".(defined $config->{mca} && $config->{mca}->get_schema_version ? $config->{mca}->get_schema_version : '?');
    
    # other version data
    my $version_data = get_version_data($config);
    $version_string .= "\n## $_ version ".$version_data->{$_} for keys %$version_data;
    
    # add key for extra column headers based on config
    my $extra_column_keys = join "\n",
        map {'## '.$_.' : '.$COL_DESCS{$_}}
        map {@{$_->{cols}}}
        grep {defined $config->{$_->{flag}}}
        @EXTRA_HEADERS;
    
    my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v$VERSION
## Output produced at $time
$db_string
## $version_string
## Extra column keys:
$extra_column_keys
HEAD
   
    $header .= join("", map {sprintf("## %s : %s\n", @{$_})} @{get_plugin_headers($config)});
    
    # add headers
    print $out_file_handle $header;
    print $html_file_handle $header if defined($config->{html});
    
    # add custom data defs
    if(defined($config->{custom})) {
        foreach my $custom(@{$config->{custom}}) {
            print $out_file_handle '## '.$custom->{name}." : ".$custom->{file}.' ('.$custom->{type}.")\n";
            print $html_file_handle '## '.$custom->{name}." : ".$custom->{file}.' ('.$custom->{type}.")\n" if defined($config->{html});
            
            foreach my $field(@{$custom->{fields}}) {
              print $out_file_handle '## '.$custom->{name}."_".$field." : ".$field." field from ".$custom->{file}."\n";
              print $html_file_handle '## '.$custom->{name}."_".$field." : ".$field." field from ".$custom->{file}."\n" if defined($config->{html});
            }
        }
    }
    
    # add column headers
    print $out_file_handle '#', (join "\t", @{$config->{fields}});
    print $out_file_handle "\n";
    
    if(defined($config->{html})) {
        print $html_file_handle html_table_headers($config, $config->{fields});
    }
    
    return $out_file_handle;
}

sub get_plugin_headers {
  my $config = shift;

  my @headers = ();

  for my $plugin (@{ $config->{plugins} }) {
    if (my $hdr = $plugin->get_header_info) {
      for my $key (sort keys %$hdr) {
        my $val = $hdr->{$key};

        push @headers, [$key, $val];
      }
    }
  }

  return \@headers;
}

# convert a variation feature to a line of output
sub convert_vf {
    my $config = shift;
    my $vf = shift;
    
    my $convert_method = 'convert_to_'.lc($config->{convert});
    my $method_ref   = \&$convert_method; 
    
    my $line = &$method_ref($config, $vf);
    my $handle = $config->{out_file_handle};
    
    if(scalar @$line) {
        print $handle join "\t", @$line;
        print $handle "\n";
    }
}

# converts to Ensembl format
sub convert_to_ensembl {
    my $config = shift;
    my $vf = shift;
    
    return [
        $vf->{chr} || $vf->seq_region_name,
        $vf->start,
        $vf->end,
        $vf->allele_string,
        $vf->strand,
        $vf->variation_name
    ];
}


# converts to pileup format
sub convert_to_pileup {
    my $config = shift;
    my $vf = shift;
    
    # look for imbalance in the allele string
    my %allele_lengths;
    my @alleles = split /\//, $vf->allele_string;
    
    foreach my $allele(@alleles) {
        $allele =~ s/\-//g;
        $allele_lengths{length($allele)} = 1;
    }
    
    # in/del
    if(scalar keys %allele_lengths > 1) {
        
        if($vf->allele_string =~ /\-/) {
            
            # insertion?
            if($alleles[0] eq '-') {
                shift @alleles;
            
                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '+'.$alleles[$i];
                }
            }
            
            else {
                @alleles = grep {$_ ne '-'} @alleles;
                
                for my $i(0..$#alleles) {
                    $alleles[$i] =~ s/\-//g;
                    $alleles[$i] = '-'.$alleles[$i];
                }
            }
            
            @alleles = grep {$_ ne '-' && $_ ne '+'} @alleles;
            
            return [
                $vf->{chr} || $vf->seq_region_name,
                $vf->start - 1,
                '*',
                (join "/", @alleles),
            ];
        }
        
        else {
            warn "WARNING: Unable to convert variant to pileup format on line number ", $config->{line_number} unless defined($config->{quiet});
            return [];
        }
        
    }
    
    # balanced sub
    else {
        return [
            $vf->{chr} || $vf->seq_region_name,
            $vf->start,
            shift @alleles,
            (join "/", @alleles),
        ];
    }
}

# converts to HGVS (hackily returns many lines)
sub convert_to_hgvs {
    my $config = shift;
    my $vf = shift;
    
    # ensure we have a slice
    $vf->{slice} ||= get_slice($config, $vf->{chr}, undef, 1);
    
    my $tvs = $vf->get_all_TranscriptVariations;
    
    my @return = values %{$vf->get_all_hgvs_notations()};
    
    if(defined($tvs)) {
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'c')}} @$tvs;
        push @return, map {values %{$vf->get_all_hgvs_notations($_->transcript, 'p')}} @$tvs;
    }
    
    return [join "\n", @return];
}

# prints a line of output from the hash
sub print_line {
    my $config = shift;
    my $line = shift;
    return 0 unless defined($line);
    
    my $output;
    my $html_fh = $config->{html_file_handle};
    
    # JSON
    if(ref($line) eq 'HASH' && defined($config->{json})) {
      $config->{json_obj} ||= JSON->new;
      $output = $config->{json_obj}->encode($line);
    }
    
    # normal
    elsif(ref($line) eq 'HASH') {
        my %extra = %{$line->{Extra}};
        
        # create extra field order?
        if(!defined($config->{field_order})) {
          my @extra_fields =
            map {@{$_->{cols}}}
            grep {defined $config->{$_->{flag}}}
            @EXTRA_HEADERS;
          
          $config->{field_order}->{$extra_fields[$_]} = $_ for 0..$#extra_fields;
        }
        
        $line->{Extra} = join ';',
          map { $_.'='.$line->{Extra}->{$_} }
          sort {$config->{field_order}->{$a} <=> $config->{field_order}->{$b}}
          keys %{ $line->{Extra} || {} };
        
        # if the fields have been redefined we need to search through in case
        # any of the defined fields are actually part of the Extra hash
        $output = join "\t", map {
            (defined $line->{$_} ? $line->{$_} : (defined $extra{$_} ? $extra{$_} : '-'))
        } @{$config->{fields}};
        
        if(defined($config->{html})) {
          my @tr_data;
          
          foreach my $field(@{$config->{fields}}) {
            my $value = $line->{$field} || $extra{$field} || '-';
            
            push @tr_data, $value eq '-' ? $value : linkify($config, $field, $value, $line);
          }
          
          print $html_fh Tr(map {td($_)} @tr_data);
        }
    }
    
    # gvf/vcf
    else {
        $output = $$line;
    }
    
    my $fh = $config->{out_file_handle};
    print $fh "$output\n";
    
    return 1;
}

sub summarise_stats {
    my $config = shift;
    
    # convert gene and transcript hashes to counts
    for my $type(qw(gene transcript regulatoryfeature)) {
      $config->{stats}->{$type} = scalar keys %{$config->{stats}->{$type}} if defined $config->{stats}->{$type};
    }
    
    # tot up chromosome counts
    foreach my $chr(keys %{$config->{stats}->{chr}}) {
      $config->{stats}->{chr_totals}->{$chr} += $config->{stats}->{chr}->{$chr}->{$_} for keys %{$config->{stats}->{chr}->{$chr}};
      
      my $start = 0;
      my %tmp;
      
      while($start <= $config->{stats}->{chr_lengths}->{$chr}) {
        $tmp{$start / 1e6} = $config->{stats}->{chr}->{$chr}->{$start} || 0;
        $start += 1e6;
      }
      
      $config->{stats}->{chr}->{$chr} = \%tmp;
    }
    
    # convert allele changes to Ts/Tv
    map {$config->{stats}->{ts_tv}->{$ts_tv{$_}} += $config->{stats}->{allele_changes}->{$_}} keys %{$config->{stats}->{allele_changes}} if defined($config->{stats}->{allele_changes});
    
    # flesh out protein_pos
    if(defined($config->{stats}->{protein_pos})) {
      if(defined($config->{stats}->{protein_pos}->{10})) {
        $config->{stats}->{protein_pos}->{9} += $config->{stats}->{protein_pos}->{10};
        delete $config->{stats}->{protein_pos}->{10};
      }
      $config->{stats}->{protein_pos}->{$_} ||= 0 for (0..9);
      
      my %tmp = map {$_.'0-'.($_+1).'0%' => $config->{stats}->{protein_pos}->{$_}} keys %{$config->{stats}->{protein_pos}};
      $config->{stats}->{protein_pos} = \%tmp;
    }
    
    # coding cons
    foreach my $con(qw(missense_variant synonymous_variant coding_sequence_variant stop_lost stop_gained frameshift_variant inframe_insertion inframe_deletion)) {
      $config->{stats}->{coding}->{$con} = $config->{stats}->{consequences}->{$con} if defined($config->{stats}->{consequences}->{$con});
    }
    
    # get ranks to sort
    my %cons_ranks = map { $_->{SO_term} => $_->{rank} } values %Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES;
    
    # create pie chart hashes
    my @charts = (
      {
        id => 'var_class',
        title => 'Variant classes',
        header => ['Variant class', 'Count'],
        data => $config->{stats}->{classes},
        type => 'pie',
        sort => 'value',
        height => 200,
      },
      {
        id => 'var_cons',
        title => 'Consequences (most severe)',
        header => ['Consequence type', 'Count'],
        data => $config->{stats}->{var_cons},
        type => 'pie',
        sort => \%cons_ranks,
        colours => $colour_keys{consequences},
      },
      {
        id => 'consequences',
        title => 'Consequences (all)',
        header => ['Consequence type', 'Count'],
        data => $config->{stats}->{consequences},
        type => 'pie',
        sort => \%cons_ranks,
        colours => $colour_keys{consequences},
      },
      {
        id => 'coding',
        title => 'Coding consequences',
        header => ['Consequence type', 'Count'],
        data => $config->{stats}->{coding},
        type => 'pie',
        sort => \%cons_ranks,
        colours => $colour_keys{consequences},
      }
    );
    
    foreach my $tool(qw(SIFT PolyPhen)) {
      my $lc_tool = lc($tool);
      
      push @charts, {
        id => $lc_tool,
        title => $tool.' summary',
        header => ['Prediction', 'Count'],
        data => $config->{stats}->{$tool},
        type => 'pie',
        height => 200,
        sort => 'value',
        colours => $colour_keys{$lc_tool},
      } if defined($config->{$lc_tool});
    }
    
    push @charts, {
      id => 'chr',
      title => 'Variants by chromosome',
      header => ['Chromosome','Count'],
      data => $config->{stats}->{chr_totals},
      sort => 'chr',
      type => 'bar',
      options => '{legend: {position: "none"}}',
    };
    
    foreach my $chr(sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$config->{stats}->{chr}}) {
      my $chr_id = $chr;
      $chr_id =~ s/\./\_/g;
      
      push @charts, {
        id => 'chr_'.$chr_id,
        title => 'Distribution of variants on chromosome '.$chr,
        header => ['Position (mb)', 'Count'],
        data => $config->{stats}->{chr}->{$chr},
        sort => 'chr',
        type => 'area',
        options => '{hAxis: {title: "Position (mb)", textStyle: {fontSize: 8}}, legend: {position: "none"}}',
        no_table => 1,
        no_link => 1,
      };
    }
    
    push @charts, {
      id => 'protein',
      title => 'Position in protein',
      header => ['Position in protein (percentile)','Count'],
      data => $config->{stats}->{protein_pos},
      sort => 'chr',
      type => 'bar',
      no_table => 1,
      options => '{hAxis: {title: "Position in protein (percentile)", textStyle: {fontSize: 10}}, legend: {position: "none"}}',
    };
    
    my @run_stats_rows = (
      ['VEP version (API)', $VERSION.' ('.$config->{reg}->software_version.')'],
      ['Cache/Database', ($config->{cache} ? $config->{dir} : ($config->{mca} ? $config->{mca}->dbc->dbname." on ".$config->{mca}->dbc->host : '?'))],
      ['Species', $config->{species}],
      ['Command line options', pre(join(" ", @{$config->{stats}->{options}}))],
      ['Start time', $config->{stats}->{start_time}],
      ['End time', $config->{stats}->{end_time}],
      ['Run time', $config->{stats}->{run_time}." seconds"],
      ['Input file (format)', $config->{input_file}.' ('.uc($config->{format}).')'],
      [
        'Output file',
        $config->{output_file}.
        (defined($config->{html}) ? ' '.a({href => $config->{output_file}.'.html'}, '[HTML]') : '').
        ' '.a({href => $config->{output_file}}, '[text]')
      ],
    );
    
    my @general_stats_rows = (
      ['Lines of input read', $config->{line_number}],
      ['Variants processed', $config->{stats}->{var_count}],
      ['Variants remaining after filtering', $config->{stats}->{filter_count}],
      ['Lines of output written', $config->{stats}->{out_count}],
      [
        'Novel / existing variants',
        defined($config->{stats}->{existing}) ?
        sprintf("%s (%.1f\%) / %s (%.1f\%)",
          $config->{stats}->{var_count} - $config->{stats}->{existing},
          100 * (($config->{stats}->{var_count} - $config->{stats}->{existing}) / $config->{stats}->{var_count}),
          $config->{stats}->{existing},
          100 * ($config->{stats}->{existing} / $config->{stats}->{var_count}),
        )
        : '-'
      ],
      ['Overlapped genes', $config->{stats}->{gene}],
      ['Overlapped transcripts', $config->{stats}->{transcript}],
      ['Overlapped regulatory features', $config->{stats}->{regulatoryfeature} || '-'],
    );
    
    # get file handle
    my $fh = $config->{stats_file_handle};
    
    # text output
    if(defined($config->{stats_text})) {
      print $fh "[VEP run statistics]\n";
      print $fh join("\t", map {s/\<.+?\>//g; $_} @{$_})."\n" for @run_stats_rows;
      
      print $fh "\n[General statistics]\n";
      print $fh join("\t", map {s/\<.+?\>//g; $_} @{$_})."\n" for @general_stats_rows;
      
      foreach my $chart(@charts) {
        print $fh "\n[".$chart->{title}."]\n";
        print $fh join("\t", ($_, $chart->{data}->{$_}))."\n" for @{sort_keys($chart->{data}, $chart->{sort})};
      }
    }
    
    # html output
    else {
      print $fh stats_html_head($config, \@charts);
      
      # create menu
      print $fh div(
        {class => 'sidemenu'},
        div(
          {class => 'sidemenu_head'},
          "Links"
        ),
        div(
          {class => 'sidemenu_body'},
          ul(
            li([
              a({href => '#masthead'}, "Top of page"),
              a({href => '#run_stats'}, "VEP run statistics"),
              a({href => '#gen_stats'}, "General statistics"),
              map {
                a({href => '#'.$_->{id}}, $_->{title})
              } grep { !$_->{no_link} } @charts,
            ])
          ),
        )
      );
      
      print $fh "<div class='main_content'>";
      
      print $fh h3({id => 'run_stats'}, "VEP run statistics");
    
      print $fh table({class => 'stats_table'}, Tr([map {td($_)} @run_stats_rows]));
      
      # vars in/out stats
      print $fh h3({id => 'gen_stats'}, "General statistics");
      print $fh table({class => 'stats_table'}, Tr([map {td($_)} @general_stats_rows]));
      
      foreach my $chart(@charts) {
        my $height = $chart->{height} || ($chart->{type} eq 'pie' ? '400' : '200');
        
        print $fh hr();
        print $fh h3({id => $chart->{id}}, $chart->{title});
        print $fh div({id => $chart->{id}."_".$chart->{type}, style => 'width: 800px; height: '.$height.'px'}, '&nbsp;');
        print $fh div({id => $chart->{id}."_table", style => 'width: 800px; height: 200px'}, '&nbsp;') unless $chart->{no_table};
      }
      
      print $fh '</div>';
      print $fh stats_html_tail();
    }
    
    $config->{stats_file_handle}->close;
}

sub stats_html_head {
    my $config = shift;
    my $charts = shift;
    
    my ($js);
    foreach my $chart(@$charts) {
      my @keys = @{sort_keys($chart->{data}, $chart->{sort})};
      
      my $type = ucfirst($chart->{type});
      
      # add colour
      if(defined($chart->{colours})) {
        my $co = 'slices: ['.join(", ", map { $chart->{colours}->{$_} ? '{color: "'.$chart->{colours}->{$_}.'"}' : '{}' } @keys).']';
        
        if(defined($chart->{options})) {
          $chart->{options} =~ s/}$/, $co}/;
        }
        else {
          $chart->{options} = "{$co}";
        }
      }
      
      # code to draw chart
      $js .= sprintf(
        "var %s = draw$type('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]), %s);\n",
        $chart->{id}.'_'.$chart->{type},
        $chart->{id}.'_'.$chart->{type},
        $chart->{title},
        $chart->{header}->[0], $chart->{header}->[1],
        join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys),
        $chart->{options} || 'null',
      );
      
      unless($chart->{no_table}) {
        
        # code to draw table
        $js .= sprintf(
          "var %s = drawTable('%s', '%s', google.visualization.arrayToDataTable([['%s','%s'],%s]));\n",
          $chart->{id}.'_table',
          $chart->{id}.'_table',
          $chart->{title},
          $chart->{header}->[0], $chart->{header}->[1],
          join(",", map {"['".$_."',".$chart->{data}->{$_}."]"} @keys)
        );
        
        # interaction between table/chart
        $js .= sprintf(
          qq{
            google.visualization.events.addListener(%s, 'select', function() {
              %s.setSelection(%s.getSelection());
            });
            google.visualization.events.addListener(%s, 'select', function() {
              %s.setSelection(%s.getSelection());
            });
          },
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
          $chart->{id}.'_'.$chart->{type},
          $chart->{id}.'_table',
        );
      }
    }
    
    my $html =<<SHTML;
<html>
<head>
  <title>VEP summary</title>
  <script type="text/javascript" src="http://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load('visualization', '1', {packages: ['corechart','table']});
  </script>
  <script type="text/javascript">
    
    function init() {
      // charts
      $js
    }
    
    function drawPie(id, title, data, options) {    
      var pie = new google.visualization.PieChart(document.getElementById(id));
      pie.draw(data, options);
      return pie;
    }
    function drawBar(id, title, data, options) {
      var bar = new google.visualization.ColumnChart(document.getElementById(id));
      bar.draw(data, options);
      return bar;
    }
    function drawTable(id, title, data) {
      var table = new google.visualization.Table(document.getElementById(id));
      table.draw(data, null);
      return table;
    }
    function drawLine(id, title, data, options) {
      var line = new google.visualization.LineChart(document.getElementById(id));
      line.draw(data, options);
      return line;
    }
    function drawArea(id, title, data, options) {
      var area = new google.visualization.AreaChart(document.getElementById(id));
      area.draw(data, options);
      return area;
    }
    google.setOnLoadCallback(init);
  </script>
  
  
  <style type="text/css">
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }
    
    a {color: #36b;}
    a.visited {color: #006;}
    
    .stats_table {
      margin: 5px;
    }
    
    tr:nth-child(odd) {
      background-color: rgb(238, 238, 238);
    }
    
    td {
      padding: 5px;
    }
    
    td:nth-child(odd) {
      font-weight: bold;
    }
    
    h3 {
      color: #666;
    }
    
    .masthead {
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }
    
    .main {
      padding: 10px;
    }
    
    .gradient {
      background: #333366; /* Old browsers */
      background: -moz-linear-gradient(left,  #333366 0%, #ffffff 100%); /* FF3.6+ */
      background: -webkit-gradient(linear, left top, right top, color-stop(0%,#333366), color-stop(100%,#ffffff)); /* Chrome,Safari4+ */
      background: -webkit-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Chrome10+,Safari5.1+ */
      background: -o-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Opera 11.10+ */
      background: -ms-linear-gradient(left,  #333366 0%,#ffffff 100%); /* IE10+ */
      background: linear-gradient(to right,  #333366 0%,#ffffff 100%); /* W3C */
      filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#333366', endColorstr='#ffffff',GradientType=1 ); /* IE6-9 */
      
      padding: 0px;
      height: 80px;
      width: 500px;
      float: right;
      display: inline;
    }
    
    .main_content {
      margin-left: 300px;
    }
    
    .sidemenu {
      width: 260px;
      position: fixed;
      border-style: solid;
      border-width: 2px;
      border-color: rgb(51, 51, 102);
    }
    
    .sidemenu_head {
      width: 250px;
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      padding: 5px;
    }
    
    .sidemenu_body {
      width: 250px;
      padding: 5px;
    }
  </style>
</head>
<body>
<div id="masthead" class="masthead">
  <div style="float: left; display: inline; padding: 10px; height: 80px;">
    <a href="http://www.ensembl.org/"><img src="http://static.ensembl.org/i/e-ensembl.png"></a>
  </div>
  
  <div style="float: right; display: inline; height: 80px; background: white; padding: 10px;">
    <a href="http://www.ensembl.org/info/docs/variation/vep/vep_script.html"><img src="http://www.ensembl.org/img/vep_logo.png"></a>
  </div>
  <div class="gradient">
  </div>
</div>
<div class="main">
SHTML

    return $html;
}

sub sort_keys {
  my $data = shift;
  my $sort = shift;
  
  my @keys;
  
  # sort data
  if(defined($sort)) {
    if($sort eq 'chr') {
      @keys = sort {($a !~ /^\d+$/ || $b !~ /^\d+/) ? $a cmp $b : $a <=> $b} keys %{$data};
    }
    elsif($sort eq 'value') {
      @keys = sort {$data->{$a} <=> $data->{$b}} keys %{$data};
    }
    elsif(ref($sort) eq 'HASH') {
      @keys = sort {$sort->{$a} <=> $sort->{$b}} keys %{$data};
    }
  }
  else {
    @keys = keys %{$data};
  }
  
  return \@keys;
}

sub stats_html_tail {
  return "\n</div></body>\n</html>\n";
}

sub html_head {
    my $txt_file = $config->{output_file};
    my $stats_file = $config->{stats_file};
    my $html =<<HTML;
<html>
<head>
  <title>VEP output</title>
  <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.8.2/jquery.min.js"></script>
  <script src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>
  <script type="text/javascript" language="javascript">
    \$(document).ready(function() {
      \$('#data').dataTable({
        "sPaginationType": "full_numbers"
      });
    });
    
    function fnShowHide( iCol ) {
      /* Get the DataTables object again - this is not a recreation, just a get of the object */
      var oTable = \$('#data').dataTable();
       
      var bVis = oTable.fnSettings().aoColumns[iCol].bVisible;
      oTable.fnSetColumnVis( iCol, bVis ? false : true );
    }
    
    function showAllCols() {
      var oTable = \$('#data').dataTable();
      for (var i=0;i<oTable.fnSettings().aoColumns.length;i++) { 
        oTable.fnSetColumnVis(i, true);
      }
    }
    
    function showHide(lyr) {
      var lyrobj = document.getElementById(lyr);
      
      if(lyrobj.style.height == "0px") {
        lyrobj.style.height = "";
        lyrobj.style.display = "";
      }
      
      else {
        lyrobj.style.height = "0px";
        lyrobj.style.display = "none";
      }
    }
  </script>
  <style type="text/css">
    \@import "http://www.datatables.net/release-datatables/media/css/demo_table.css";
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }
    
    a {color: #36b;}
    a.visited {color: #006;}
    
    th {
      font-size: 11px;
    }
    td {
      font-size: 11px;
    }
    
    .masthead {
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }
    
    .main {
      padding: 10px;
    }
    
    .gradient {
      background: #333366; /* Old browsers */
      background: -moz-linear-gradient(left,  #333366 0%, #ffffff 100%); /* FF3.6+ */
      background: -webkit-gradient(linear, left top, right top, color-stop(0%,#333366), color-stop(100%,#ffffff)); /* Chrome,Safari4+ */
      background: -webkit-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Chrome10+,Safari5.1+ */
      background: -o-linear-gradient(left,  #333366 0%,#ffffff 100%); /* Opera 11.10+ */
      background: -ms-linear-gradient(left,  #333366 0%,#ffffff 100%); /* IE10+ */
      background: linear-gradient(to right,  #333366 0%,#ffffff 100%); /* W3C */
      filter: progid:DXImageTransform.Microsoft.gradient( startColorstr='#333366', endColorstr='#ffffff',GradientType=1 ); /* IE6-9 */
      
      padding: 0px;
      height: 80px;
      width: 500px;
      float: right;
      display: inline;
    }
  </style>
  </head>
  <body>
  <div id="masthead" class="masthead">
    <div style="float: left; display: inline; padding: 10px; height: 80px;">
      <a href="http://www.ensembl.org/"><img src="http://static.ensembl.org/i/e-ensembl.png"></a>
    </div>
    
    <div style="float: right; display: inline; height: 80px; background: white; padding: 10px;">
      <a href="http://www.ensembl.org/info/docs/variation/vep/vep_script.html"><img src="http://www.ensembl.org/img/vep_logo.png"></a>
    </div>
    <div class="gradient">
    </div>
  </div>
  <div class="main">
  <p>
    View: <a href="$stats_file">Summary statistics</a> |
    <a href="$txt_file">as text</a> |
    <a href="javascript:void();" onclick="showHide('header')">Show/hide header</a> |
    <a href="javascript:void();" onclick="showAllCols()">Restore columns</a>
  </p>
  <hr/>
  <pre id="header" style="height:0px; display: none;">
HTML
  return $html;
}

sub html_table_headers {
  my $config = shift;
  my $cols = shift;
  
  my @cols_copy = @$cols;
  
  my $html = qq{</pre><table id="data" class="display"><thead><tr>};
  
  $config->{_th} = join("", map {
    $cols_copy[$_] =~ s/\_/ /g;
    '<th>'.$cols_copy[$_].' '.a({href => 'javascript:void();', onclick => "fnShowHide($_);"}, img({src => 'http://www.ensembl.org/i/16/cross.png', height => 6, width => 6, style => 'border: 1px solid gray; padding: 1px;'})).'</th>'
  } (0..$#cols_copy));
  
  $html .= $config->{_th};
  $html .= qq{</thead></tr><tbody>};
  
  return $html;
}

sub linkify {
  my $config = shift;
  my $field  = shift;
  my $string = shift;
  my $line   = shift;
  
  my $species = ucfirst($config->{species});
  
  # Ensembl genes
  $string =~ s/(ENS.{0,3}G\d+|CCDS\d+\.?\d+?|N[MP]_\d+\.?\d+?)/a({href => "http:\/\/www.ensembl.org\/$species\/Gene\/Summary\?g=$1", target => "_blank"}, $1)/ge;
  
  # Ensembl transcripts
  $string =~ s/(ENS.{0,3}T\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Transcript\/Summary\?t=$1", target => "_blank"}, $1)/ge;
  
  # Ensembl regfeats
  $string =~ s/(ENS.{0,3}R\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Regulation\/Summary\?rf=$1", target => "_blank"}, $1)/ge;
  
  # variant identifiers
  $string =~ s/(rs\d+|COSM\d+|C[DMIX]\d+)/a({href => "http:\/\/www.ensembl.org\/$species\/Variation\/Summary\?v=$1", target => "_blank"}, $1)/gie;
  
  # split strings
  $string =~ s/([,;])/$1 /g;
  
  # locations
  while($string =~ m/(^[A-Z\_\d]+?:[1-9]\d+)(\-\d+)?/g) {
    my $loc = $1.($2 ? $2 : '');
    my ($chr, $start, $end) = split /\-|\:/, $loc;
    $end ||= $start;
    
    # adjust +/- 1kb
    my $view_start = $start - 10;
    my $view_end   = $end + 10;
    my $allele = $line->{Allele} || 'N';
    
    my $url =
      "http://www.ensembl.org/$species/Location/View?".
      "r=$chr:$view_start\-$view_end;format=vep_input;".
      "custom_feature=$chr%20$start%20$end%20$allele%201;".
      "custom_feature=normal";
    
    my $link = a({href => $url, target => "_blank"}, $string);
    $string =~ s/$loc/$link/;
  }
  
  return $string;
}

# outputs usage message
sub usage {
    my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version $VERSION
by Will McLaren (wm2\@ebi.ac.uk)

Help: dev\@ensembl.org , helpdesk\@ensembl.org
Twitter: \@ensembl , \@EnsemblWill

http://www.ensembl.org/info/docs/tools/vep/script/index.html

Usage:
perl variant_effect_predictor.pl [--cache|--offline|--database] [arguments]

Basic options
=============

--help                 Display this message and quit

-i | --input_file      Input file
-o | --output_file     Output file
--force_overwrite      Force overwriting of output file
--species [species]    Species to use [default: "human"]
                       
--everything           Shortcut switch to turn on commonly used options. See web
                       documentation for details [default: off]                       
--fork [num_forks]     Use forking to improve script runtime

For full option documentation see:
http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html

END

    print $usage;
}
