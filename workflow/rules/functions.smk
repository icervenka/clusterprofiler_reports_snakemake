def get_file(wildcards):
    file = Metadata[Metadata.contrast == wildcards.contrast].file
    file = INPUT_DIR + file
    return file

def get_contrast(wildcards):
    c = Metadata[Metadata.contrast == wildcards.contrast].contrast
    return c

def get_sample_expression(wildcards):
    c = Metadata[Metadata.contrast == wildcards.contrast].sample_expression
    return INPUT_DIR + c

def get_sample_sheet(wildcards):
    c = Metadata[Metadata.contrast == wildcards.contrast].sample_sheet
    return INPUT_DIR + c

def get_group_order(wildcards):
    c = Metadata[Metadata.contrast == wildcards.contrast].group_order
    return c

# def get_meta_column(column):
#     lambda wildcards: Metadata[Metadata.contrast == wildcards.contrast][column]

# #c = Metadata[Metadata.contrast == wildcards.contrast].sample_sheet
# get_contrast = get_meta_column("contrast")
# get_sample_expression = get_meta_column("sample_expression")
# get_sample_sheet = get_meta_column("sample_sheet")
# get_group_order = get_meta_column("group_order")