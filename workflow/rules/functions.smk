def get_file(wildcards):
    file = Metadata[Metadata.contrast == wildcards.contrast].file
    file = INPUT_DIR + file
    return file

def get_contrast(wildcards):
    c = Metadata[Metadata.contrast == wildcards.contrast].contrast
    return c

def create_directory(name):
    try:
        os.mkdir(OUTPUT_DIR)
    except OSError as e:
        pass