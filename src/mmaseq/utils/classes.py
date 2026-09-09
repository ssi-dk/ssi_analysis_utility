from dataclasses import dataclass, field
from pathlib import Path
from yaml import safe_load
import pandas as pd
from itertools import product

# Classes
@dataclass
class Input:
    read1: Path | None = None
    read2: Path | None = None
    assembly: Path | None = None
    config: str = "default.yaml"

    def __post_init__(self):
        if not self.read1:
            self.read1 = None
        else:
            self.read1 = Path(self.read1)

        if not self.read2:
            self.read2 = None
        else:
            self.read2 = Path(self.read2)

        if pd.isnull(self.assembly):
            self.assembly = None
        else:
            self.assembly = Path(self.assembly)

    @property
    def read_type(self) -> str | None:
        # Define paired read_type
        if self.read1 is not None and self.read2 is not None:
            return "paired"

        # Define long read type
        if self.read1 is not None:
            return "long"
        
        # Define missing reads
        return None


@dataclass
class Module:
    options: str = ""
    database: list[str] = field(default_factory=list)
    assembler: list[str] = field(default_factory=list)
    ignore: bool = False
    reads: bool = False
    config: dict = field(default_factory=dict)
    raw_results: Path | None = None
    results_final: Path | None = None


@dataclass
class Sample:
    name: str
    config_file: Path
    inputs:  Input
    outpath: Path
    modules: dict[str, Module] = field(default_factory = dict)


    def module(self, name: str) -> Module | None:
        return self.modules.get(name)


    # def show_read1(self):
    #     return self.inputs.read1
    

    # def show_read2(self):
    #     return self.inputs.read2

    
    # def show_assembly(self):
    #     return self.inputs.assembly


    # def show_module(self, module):
    #     return self.modules.get(module)


    # def show_options(self, module):
    #     return self.show_module(module).options


    # def show_database(self, module):
    #     return self.show_module(module).database


    # def show_assembler(self, module):
    #     return self.show_module(module).assembler


    # def show_config(self, module):
    #     return self.show_module(module).config


    # def show_raw_results(self, module):
    #     return self.show_module(module).raw_results


    # def show_results_final(self, module):
    #     return self.show_module(module).results_final


    def populate_modules(self):
        modules = {}

        static_fields = [
            "options",
            "database",
            "assembler",
            "ignore",
            "reads"
        ]

        with open(self.config_file) as cfg_read:
            configs_raw = safe_load(cfg_read)

        for module, components in configs_raw.items():
            opt = check_config_fields(components = components, field = "options", default = "")
            db = as_list(check_config_fields(components = components, field = "database", default = []))
            asm = as_list(check_config_fields(components = components, field = "assembler", default = []))
            ignr = check_config_fields(components = components, field = "ignore", default = False)
            rds = check_config_fields(components = components, field = "reads", default = False)

            config = {
                key: value
                for key,value in components.items()
                if key not in static_fields
            }

            # Determine output file exts
            exts = determine_outfile_exts(database = db, assembly = asm)

            # Determine sample read type
            read_type = self.inputs.read_type

            read_str = ""
            if rds and read_type:
                read_str = f"{read_type}/"
            
            raw_results = None
            results_final = None
            if not ignr:
                raw_results = [f"{self.outpath}/{read_str}{module}/{module}{ext}" for ext in exts]
                results_final = [f"{raw_res.replace('/raw', '')}" for raw_res in raw_results]


            mod = Module(
                options = opt,
                database = db,
                assembler = asm,
                ignore = ignr,
                reads = rds,
                config = config,
                raw_results = raw_results,
                results_final = results_final
            )

            modules.update({module : mod})

        self.modules = modules

        return self


## Initiation scripts

def check_config_fields(components, field, default):
    value = default
    if field in components:
        value = components.get(field)

    return value


def as_list(value):
    if not value:
        return []

    if isinstance(value, list):
        return value

    return [value]


def determine_outfile_exts(database, assembly):

    db_strings = [f"_{db}" for db in database] if database else []
    asm_strings = [f"_{asm}" for asm in assembly] if assembly else []

    # Neither database nor assembler specified
    if not db_strings and not asm_strings:
        return [".tsv"]

    # Combine every database with every assembler
    if db_strings and asm_strings:
        return [
            f"{db}{asm}.tsv"
            for db, asm in product(db_strings, asm_strings)
        ]

    # Only databases
    if db_strings:
        return [f"{db}.tsv" for db in db_strings]

    # Only assemblers
    return [f"{asm}.tsv" for asm in asm_strings]


def import_dataset(samplesheet, config_dir, outdir):
    samples = {}
    for sample in samplesheet.index:
        read1 = samplesheet.at[sample, "read1"]
        read2 = samplesheet.at[sample, "read2"]
        assembly = samplesheet.at[sample, "assembly"]

        config = Path(samplesheet.at[sample, "config"])

        config_file = config_dir / config

        if not config_file.exists():
            raise FileNotFoundError(f"Configuration file not found - Check your samplesheet: {config_file}.")

        smpl = Sample(name = sample, config_file = config_file, inputs = Input(read1, read2, assembly), outpath = outdir / sample / "raw")

        smpl.populate_modules()
        samples.update({sample: smpl})


    return samples


# # Testing    

# samplesheet = pd.read_csv(
#     "/home/cucumbergebt/micromamba/envs/MMAseq/lib/python3.14/site-packages/mmaseq/Deploy/MMAseq_Test/samplesheet_test_resolved.tsv",
#     sep = "\t",
#     na_filter=False
# ).set_index(
#     "sample_name"
# )

# config_dir = Path("/home/cucumbergebt/micromamba/envs/MMAseq/lib/python3.14/site-packages/mmaseq/config/species_configs")

# outdir = Path("/home/cucumbergebt/micromamba/envs/MMAseq/lib/python3.14/site-packages/mmaseq/Deploy/MMAseq_Test")

# samples = import_dataset(samplesheet, config_dir, outdir)

# sample = samples.get("Ec_Test")
# print(f"Sample {sample.name} with read types {sample.inputs.read_type} with {sample.inputs.read1} and maybe {sample.inputs.read2}")

# sample.module("amrfinder").raw_results
