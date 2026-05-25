import alphaquant.config.variables as aqvariables


class TestSetInputConfig:
    def setup_method(self):
        self._orig_type = aqvariables.INPUT_TYPE
        self._orig_dict = aqvariables.CONFIG_DICT

    def teardown_method(self):
        aqvariables.INPUT_TYPE = self._orig_type
        aqvariables.CONFIG_DICT = self._orig_dict

    def test_sets_globals(self):
        aqvariables.set_input_config("diann_fragion", {"key": "value"})
        assert aqvariables.INPUT_TYPE == "diann_fragion"
        assert aqvariables.CONFIG_DICT == {"key": "value"}

    def test_overwrite(self):
        aqvariables.set_input_config("type_a", {"a": 1})
        aqvariables.set_input_config("type_b", {"b": 2})
        assert aqvariables.INPUT_TYPE == "type_b"
        assert aqvariables.CONFIG_DICT == {"b": 2}
