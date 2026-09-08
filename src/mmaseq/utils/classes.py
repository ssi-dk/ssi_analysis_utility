from dataclasses import dataclass, field
from pathlib import Path

# Content classes

@dataclass
class Input:
    read1: Path | None = None
    read2: Path | None = None
    assembly: Path | None = None
    config: str

    def __post_init__(self):
        if self.read1 is not None:
            self.read1 = Path(self.read1)

        if self.read2 is not None:
            self.read2 = Path(self.read2)

        if self.assembly is not None:
            self.assembly = Path(self.assembly)

    @property
    def read_type(self) -> str | None:
        # Define paired read_type
        if self.read1 is not None and self.read2 is not None:
            return "paired"

        # Define long read type
        elif self.read1 is not None and self.read2 is None:
            return "long"
        
        # Define missing reads
        return None


    def config_file(self) -> str | None:
        ASK GPT

CONTINUE GPT

@dataclass
class Module:
    name: str
    options: str = ""
    databases: list[str] = field(default_factory=list)
    assemblers: list[str] = field(default_factory=list)
    ignore: bool = False
    reads: bool = False
    config: dict = field(default_factory=dict)


@dataclass
class Sample:
    name: str
    input:  Input
    modules: dict[str, Module] = field(default_factory = dict)


sample = mmsample(
    sample_name = "test_123",
    input = sample_input(
        read1 = "sample1_R1.fastq.gz",
        read2 = "sample2_R2.fastq.gz"
    ),
    modules = {
        "Resfinder": "this is a test",
        "serovar_detector": {"This": "is a nested dick!"}
    }
)