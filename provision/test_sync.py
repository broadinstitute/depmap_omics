from dump_workspace import dumb_wdl_parse
from sync import MethodConfig, _sync_config, WDLCache



def test_sync_config():
    wdl = """
version 1.0

task hello {
  input {
    String name
  }

  command {
    echo 'hello ${name}!'
  }
  output {
    File response = stdout()
  }
  runtime {
   docker: 'ubuntu:latest'
  }
}

workflow test {
  call hello
}
    """
    config =  MethodConfig(name="test_sync_config", inputs={}, outputs={}, wdl=dumb_wdl_parse(wdl), imports={}, wdl_sha256="x131", entity_root_type=None, prerequisites={})
    cache = WDLCache()
    _sync_config("broad-institute-achilles", "pmontgom-test", "pmontgom-test", config, cache)