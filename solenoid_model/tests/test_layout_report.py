from solenoid_model.report import generate_layout


def test_generate_section_writes_a_png():
    path = generate_layout.generate_section()
    assert path.exists()
    assert path.suffix == ".png"
    assert path.stat().st_size > 1000


def test_generate_architecture_writes_a_png():
    path = generate_layout.generate_architecture()
    assert path.exists()
    assert path.suffix == ".png"
    assert path.stat().st_size > 1000


def test_importing_does_not_mutate_global_rcparams():
    """中文字型必須以 rc_context 區域套用，不得污染全域。"""
    import matplotlib
    before = list(matplotlib.rcParams["font.sans-serif"])
    import importlib
    importlib.reload(generate_layout)
    assert list(matplotlib.rcParams["font.sans-serif"]) == before
