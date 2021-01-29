I have to say that when I had planned this article, I was not expecting the craziness of the next 6 months.
I have now spent 1.5 years at the Broad. But 6 of them have actually been remote as I worked from France. I now expect to spend the entirety of this second year, remote.
I am thus going to talk about the first non-remote, pre-Covid year.

## Finishing onboarding.
I was warned that the onboarding might last a couple of months. Being a foreigner did not help and I can finally say that now, a year later, I feel fully onboarded to the Broad and the Cambridge way of life.
Cambridge way of life often implies, amongst other things: grocery shopping in Union Square Market, getting a couple of doughnuts from Union Square Doughnuts, drinking fancy beers at Remnant or Lamplighter, asking people if they are doing a postdoc or a PhD the first time you meet, getting your baguette at Tate and enjoying the many ice creams spots Cambridge has to offer
What did it take me so long to feel perfectly unboarded? Well, except from getting up to speed with the lab's specific softwares and technologies, getting a sense of the myriads of projects  and collaborations happening around you, it's getting to know the history of projects that came before you, the reasons for things. It's getting to know power dynamics!  It is so incredibly important to know the power dynamics. I really feel that once you do, you know you are truly on board. You know you have been onboarded when you know what actions you can do and what actions you can't. Then you are able to navigate your work with much more ease. And last but not least, once you get your routine: waiting for the CC&E meetings (cell circuitry and epigenomics), the [MIA talks](https://www.youtube.com/watch?v=MKD1ScU_XJs&list=PLlMMtlgw6qNjROoMNTBQjAcdx53kV50cS) (models, inference and algorithms) .
Once you are onboarded this is also when you start getting ideas for projects you could launch and whom to collaborate with.

## What I have done
During this year I have been able to present papers and lead team meetings.
I have been able to present my research during the yearly Broad research conference, an event with 3000 people, covering all imaginable research in molecular biology .
I was also able to interact with many teams across the Broad for small projects, either helping them solve issues or getting help from them to solve my issues.
As you may know from other posts: _sesee [AML project](https://jkobject.com/projects/predicting-dependencies-from-enhancers/), [depmap omics](https://jkobject.com/projects/depmap-omics-howto/) I am the principal bio-informaticist of two main projects: one at the Broad and one at the Dana-Farber Cancer Institute.
However, I was also able to participate in a few side projects as well (side projects for me :wink:) regarding specific cancer dependencies.

### Diversity of Projects
One of those projects was focused on the WRN dependency in MSI cancers, see [wrn dependencies](https://jkobject.com/projects/wrn-dependencies/), and one was on the VPS4A dependency in VPS4B deficient cancers, see [vps4a dependency](). These projects were driven primarily by biologists and experimentalists. They came to us because they needed some additional computational analysis.
Once the goal of the project was understood, this was "fairly trivial" analysis, wrapped up in less than 2 weeks. What takes 90% of my time are the two main projects.
For the first 8 months of my time I was also involved in revamping [CCLF](https://cellfactory.broadinstitute.org/)'s computational pipeline. I also laid down a first SOP and road-map for the data processing and release .

### CCLE
For CCLE, a set of circumstances made it so that I ended up being the only computationalist on the project during a period of 10 months. Hopefully Neekesh Dharia, my mentor at the time was able to provide the biology intuition in all of that. I went on to maintain, QC and productionalize the pipeline. My main focus was on reproducibility. It quickly appeared that nothing of what we were producing could be reproduced by anyone. Even ourselves. We had to work with many people from [GtEX](https://gtexportal.org/home/), [CGA](https://www.broadinstitute.org/cancer/cancer-genome-computational-analysis), [DSP](https://www.broadinstitute.org/data-sciences-platform), [GP](https://www.broadinstitute.org/reading-and-editing-biology/genomics-platform).. At the end, I layed down a plan taking from everything we had learned along the way.

#### CCLE: Main Issues
The two major issues faced in our work have been on the data sharing process. Even cell lines, released at this scale, need to have a regulatory framework built around them. Issues are mostly around pseudo-somatic[1] and [germline](https://en.wikipedia.org/wiki/Germline_mutation) mutations and which ones can be released to whom, and how?
The second one is mainly regarding the size of the project. This project is between a large science experiment and a fully fledged product. The amount of users that we have and the sensitivity of their projects make any of our changes have big consequences downstream for a lot of human cancer research\*. 

This entails that any change we want to make needs to be carefully reviewed and thought through. This makes for a lot of frustrations. Things can be slow. Mistakes & non optimal decisions can live on for much longer because the effort needed to change them becomes very high. 

\* _I think that, at DepMap, we all handle this pressure fairly well. This drives us to put forward the best of ourselves while still being ok with our mistakes. We know that everyone makes mistakes and that the state of genetic cancer research and computational biology as a whole does not reach anyone's expectations, yet._

#### CCLE: What is to come
Most of the main objectives are still works in progress. We are in the process of hiring a new ACB to take on most of the workload around improving, productionalizing and running the pipeline itself. :construction:
More recently, we have switched to WGS -whole genome- data which required us to rethink the flow data architecture.
My main focus  for the coming month will be to have it in a good state of documentation and to have clearly laid down SOPs and directions of improvements for the new ACB -entry level computational biologist-.
I will also soon start to focus full time on a machine learning project (undisclosed for now), related to CCLE and Dependencies.

### AMLproject
Here we are working mainly with very deep sequencing of a handful of cell lines under many different conditions . The dynamics are very different:
- we are pursuing a very specific goal  in understanding a biological phenomena
- we are all young scientist learning as we go and doing more cutting edge experiments (and thus, analysis)
- reliability and efficiency is really a side goal (until it can't be)
- many, many different types of experiments of various quality exist.
We have now almost proved what we need to in order to publish this work in a high tier journal. Maxim Pimkin has really been the master mind  behind everything related to the biology in this project. His intuition  has led most of our decisions, which is also a big difference compared to CCLE. Other very talented people have worked hard on this project, but I like to think of myself as the lead computational biologist.

#### AML: Main Issues
Something that is notoriously hard in a project is to work with low foresight. Here each new experiment would trigger new unforeseen ones. 
In code there is two ways to do things: 
- iteratively: 1st version, 2nd version etc.. --> going toward a metric of quality
- in one go: you have to build something and you build it block by block
Here we proceeded with the first method, however we can never know what the next version will be and we could go as much toward an improved version than toward a version with new functionalities.
This obviously takes time. Each version is also very specific toward one specific thing. Transmitting this information and approach is particularly hard for biologists. They often expect the latter. It is also hard to estimate the time things will take. A very similar experiment could run much faster than a slightly different one for example.
A marginal amount of work on a first experiment could save countless hours in the coming ones, etc..
A very similar experiment could run much faster than a slightly different one for example.

Interacting across disciplines is hard but leads to some of the best projects, and we can see some examples here. :sparkles:

#### AML: What is to come
During this project we have discovered many interesting properties of the leukemia regulatory circuitry. We have also discovered interesting features of [slam-seq](https://www.lexogen.com/slamseq-metabolic-rna-labeling/) . We want it to become an additional side -letter format- paper. Moreover we have left open many biological questions.
One of our main objectives is quite similar to CCLE's. We need an additional ACB who will learn from me and learn how to use my pipelines. This will also free me from wrapping up some analysis and thinking more about the paper and how to better release the code. 
But I want to also try and be able to finish our analysis of large scale [transcription factor](https://en.wikipedia.org/wiki/Transcription_factor) co-binding. 
Let's meet in a month to see how things have turned out! :wave:

---
[1] I call pseudo-somatic, [somatic](https://en.wikipedia.org/wiki/Somatic_mutation) mutations in cell lines as they are somewhat more losely defined than somatic in cancer tissues, they regroup immortalizing mutations -due to the process of [immortalizing a cell line](https://en.wikipedia.org/wiki/Immortalised_cell_line), [private](https://www.medicinenet.com/script/main/art.asp?articlekey=5048) snps and cancer mutations that also exist as germlines.