def pytest_addoption(parser):
    parser.addoption("-P", "--path", action="store", default="./")
    parser.addoption("-L", "--launcher_cmd", action="store", default="")
    parser.addoption("--launcher_opts", action="store", default="")
    parser.addoption("--node_range", action="store", default="")
    parser.addoption("-M", "--implementation_mask", action="store", default="31")

def pytest_generate_tests(metafunc):
    option_value = metafunc.config.option.path
    if 'path' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("path",[option_value])

    option_value = metafunc.config.option.launcher_cmd
    if 'launcher_cmd' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("launcher_cmd",[option_value])
        
    option_value = metafunc.config.option.launcher_opts
    if 'launcher_opts' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("launcher_opts",[option_value])

    option_value = metafunc.config.option.node_range
    if 'node_range' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("node_range",[option_value])

    option_value = metafunc.config.option.implementation_mask
    if 'implementation_mask' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("implementation_mask",[option_value])
