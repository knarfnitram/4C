baci_auto_define_module()

set(_dependencies
    # cmake-format: sortable
    adapter
    core
    fluid
    )
baci_add_dependency(${AUTO_DEFINED_MODULE_NAME} ${_dependencies})
