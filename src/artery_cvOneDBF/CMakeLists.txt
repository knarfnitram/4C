baci_auto_define_module()

set(_artery_cvOneDBF
    # cmake-format: sortable
    adapter
    core
    fluid
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_artery_cvOneDBF})
