from .utils import parse_setup

args = parse_setup()

    # if test:
    #     logger.info("Test run initiated. Will ignore irrelevant user arguments!")
    #     config = f"{PKG_CONFIGS}/Test.yaml"
    #     samplesheet_file = f"{DATA_DIR}/samplesheet.tsv"
    #     outdir = outdir / "Test"
    #     rules.append("all")

    #     # Fix: normalize paths in test
    #     if resolve:
    #         resolve_samplesheet_paths(samplesheet_file)

    #     config = create_config(samplesheet_file, 
    #                            outdir, 
    #                            deploy_dir, 
    #                            config, 
    #                            force = True)

    #     link_assemblies(samplesheet_file, 
    #                     SPE_CONFIGS, 
    #                     outdir)

    #     command = create_command(threads, 
    #                              config, 
    #                              conda_dir, 
    #                              arguments, 
    #                              rules)


        # if samplesheet_file is None:
        #     logger.error("Samplesheet file user argument missing. Aborting!")
        #     sys.exit(1)