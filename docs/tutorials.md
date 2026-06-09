# Tutorials 


## Deployment of the conda environments,databases and run the test set

After installing the pipeline, the next step is to deploy the conda environments and databases that are required for running the analyses. This can be done using the `mmadeploy` command, which will set up the necessary environments and databases in a specified directory. The deployment process ensures that all the required tools and databases are available and properly configured for use with the pipeline. Moreover, it will run the test set to verify that everything is correctly installed and configured. Running the test set is a good way to verify that the pipeline is correctly installed and configured on your system. It will execute a predefined set of analyses using example data, and it will help you familiarize with the expected input and output formats. 

Verify that the pipeline was correctly installed:
```bash
mmaseq --help
```

To run the test set, simply execute the following command:
```bash
mmadeploy --deploy_dir path/to/deploy_dir  --threads 4
```


## Run the pipeline on your data

To run the pipeline on your own data, you will need to prepare a samplesheet in a tab-separated format (TSV)  that contains the necessary information about your samples, following the format provided in [Samplesheet Configuration](installation.md#samplesheet-configuration). Once you have prepared your samplesheet, you can run the pipeline using the `mmaseq` command, specifying the path to your samplesheet, the deployment directory where the conda environments and databases are located, and the output directory where you want the results to be saved.

Run the pipeline with the following command:

```bash
mmaseq --samplesheet path/to/samplesheet.tsv 
       --deploy_dir path/to/deploy_dir 
       --outdir path/to/output_dir
```

