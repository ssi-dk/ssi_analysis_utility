import argparse
from .paths import *

def parse_setup():
    parser = argparse.ArgumentParser(
        description = """
        Run pipeline on inbuild test set and setup module environments and databases.
        """
    )

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = PKG_DIR / "Deploy",
        help = f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default {PKG_DIR}/Deploy)
        """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        default = CWD / "MMAseq_Results",
        help = """
            Directory used for test dataset and analysis results. (Default ./MMAseq_Results)
        """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = """
            Amount of threads (cores) to dedicate for executing the pipeline. 
            (Default 4)
        """
    )

    parser.add_argument(
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
            Add debug messages during execution. Mostly used for development 
            and debugging purposes
        """
    )

    return parser.parse_args()


def parse_create():
    parser = argparse.ArgumentParser(
        description = """
        Create samplesheet for execution
        """
    )

    parser.add_argument(
        "--indir",
        dest = "indir",
        default = None,
        help = """
            Input directory MUST be specified if the samplesheet does 
            not yet exist. (Default None)
            Input directory will be screened for `.fasta` and `fastq.gz` 
            files, sample_names will be infered from the detected files, 
            and used to populate a samplesheet. After samplesheet creation, 
            the pipeline will be executed in dry-run mode (simulated run)
        """
    )

    parser.add_argument(
    "--outdir",
    dest = "outdir",
    default = CWD / "MMAseq_Results",
    help = """
        Directory used for storing the samplesheet. (Default ./MMAseq_Results)
        Samplesheet will be stored in ourdir as 'samplesheet.tsv'
    """
    )

    parser.add_argument(
    "--debug",
    dest = "debug",
    action = "store_true",
    help = """
        Add debug messages during execution. (Default False) 
        Mostly used for development and debugging purposes.
    """
    )

    return parser.parse_args()


def parse_mmaseq():
    parser = argparse.ArgumentParser(
        description = """
        Execute Mixed Microbial Analysis on Sequencing data (MMAseq)
        """
    )

    parser.add_argument(
        "--samplesheet",
        dest = "samplesheet_file",
        required = True,
        help = """
            Path to samplesheet TSV used by the pipeline. 
            If the samplesheet doesn't exist, `indir` must be 
            specified to create one.
        """
    )
    

    parser.add_argument(
        "--deploy_dir",
        dest = "deploy_dir",
        default = PKG_DIR / "Deploy",
        help = f"""
            Directory used to deploy virtual environment and databases 
            used during pipeline execution. To reinstall environments 
            and/or databases, remove the `conda/` and/or the `Databases/` 
            folders in the deployment directory. (Default {PKG_DIR}/Deploy)
        """
    )

    parser.add_argument(
        "--outdir",
        dest = "outdir",
        default = CWD / "MMAseq_Results",
        help = """
            Directory used for storing analysis results. (Default ./MMAseq_Results)
        """
    )

    parser.add_argument(
        "--threads",
        dest = "threads",
        default = 4,
        help = """
            Amount of threads (cores) to dedicate for executing the pipeline. 
            (Default 4)
        """
    )

    parser.add_argument(
        "--config",
        dest = "config",
        default = PKG_CONFIGS / "config.yaml",
        help = f"""
            Configuration file location. (Default {PKG_CONFIGS}/config.yaml)
        """
    )

    # parser.add_argument(
    #     "--test",
    #     dest = "test",
    #     action = "store_true",
    #     help = """
    #         Perform test run of all modules. Can be used to install all 
    #         conda environments and databases. The config file will be 
    #         generated as [CONFIG_DIR]/Test.yaml and output will be stored 
    #         in a Test/ folder of the specified output directory. 
    #         (Default ./MMAseq_Results/Test)
    #     """
    # )

    parser.add_argument(
        "--resolve",
        dest = "resolve",
        action = "store_true",
        help = """
            Resolves absolute paths from samplesheet columns, will 
            overwrite samplesheet. (Default False)
        """
    )

    parser.add_argument(
        "--ignore_assemblies",
        dest = "ignore_assemblies",
        action = "store_true",
        help = """
            Avoid creating symbolic links of the assemblies into the 
            pipeline output directory. (Default False) 
            If this is not specified, symbolic links will be created 
            to the output directory, hence avoidubg assembly steps in the 
            pipeline. Use this option to enforce the pipeline to create 
            assemblies (Will take extra time) rather than relying on those 
            specified in the samplesheet.
        """
    )

    parser.add_argument(
        "--debug",
        dest = "debug",
        action = "store_true",
        help = """
            Add debug messages during execution. (Default False) 
            Mostly used for development and debugging purposes.
        """
    )

    return parser.parse_args()