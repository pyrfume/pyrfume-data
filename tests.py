import unittest

import pyrfume
import toml


class ManifestTestCase(unittest.TestCase):
    """Unit tests for the manifests"""
    
    def setUp(self):
        self.archives = pyrfume.list_archives()

    def test_load(self):
        """Test backends."""
        errors = []
        for archive in self.archives:
            try:
                manifest = pyrfume.load_manifest(archive)
            except toml.TomlDecodeError as e:
                errors.append('%s: %s' % (archive, e))
        if errors:
            self.fail('\n'.join(errors))

