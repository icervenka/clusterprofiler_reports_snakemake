import glob
# ruleorder: move_libs > remove_libs > change_libpath
#
# rule move_libs:
#     input:
#         directory(REPORT_OUTDIR + Metadata.contrast[0] + "/libs")
#     output:
#         directory(REPORT_OUTDIR + "libs")
#     shell:
#         "mv {input} {output}"
#
# rule remove_libs:
#     input:
#         expand(REPORT_OUTDIR + "{contrast}" + "/libs", contrast = Metadata.contrast)
#     output:
#         touch('libs.removed')
#     shell:
#         "rm -r {input}"
# rule change_libpath_completed:
#     input:
#         "libpath_changed"

rule change_libpath:
    input:
        #expand(REPORT_OUTDIR + "{contrast}/{type}_unique.html", contrast = Metadata.contrast, type = types)
        glob_wildcards(REPORT_OUTDIR + "{html}")
    output:
        touch("libpath_changed")
    wildcard_constraints:
        html = '\w+\.html$'
    run:
        #items = glob.glob(input[0])
        # print(input)
        # print(input[0])
        # for item in input:
        print(input[0])
        i = open(input[0])
        o = open(input[0], "w")
        o.write(i.read().replace('="libs/', '="../libs/'))

# rule move_libs_all:
#     input:
#         'libpath.changed'
