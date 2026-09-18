"""
Every shipped config must validate. This runs in CI, so a config that has
drifted out of step with the validator fails the build rather than failing
someone's run three weeks later.
"""

import os

import pytest
import yaml

from lib.config import validate, ConfigError, load


def test_at_least_one_config_is_shipped(configs):
    assert configs, "no pathogen configs found in config/pathogen/"


def test_every_shipped_config_validates(configs, repo_root):
    os.chdir(repo_root)
    failed = {}
    for p in configs:
        try:
            load(p)
        except ConfigError as e:
            failed[p.name] = str(e)
    assert not failed, f"configs failed validation: {failed}"


def test_template_refuses_to_run_unedited(repo_root):
    """
    The template must NOT validate. A placeholder taxid that passes silently
    would let someone run the template unchanged and download nothing.
    """
    os.chdir(repo_root)
    t = repo_root / "config" / "pathogen" / "_template.yaml"
    cfg = yaml.safe_load(t.read_text())
    cfg["_path"] = str(t)
    errors, _ = validate(cfg)
    assert errors, "the template validated; it must not"


def test_config_ids_are_unique(configs):
    ids = [yaml.safe_load(p.read_text())["id"] for p in configs]
    assert len(set(ids)) == len(ids), f"duplicate pathogen ids: {ids}"


def test_config_id_matches_filename(configs):
    """builds/<id>/ is derived from id; a mismatch makes builds hard to find."""
    for p in configs:
        assert yaml.safe_load(p.read_text())["id"] == p.stem


@pytest.mark.parametrize("required", ["title", "build_name", "maintainer"])
def test_nextstrain_block_is_complete(configs, required):
    for p in configs:
        ns = yaml.safe_load(p.read_text()).get("nextstrain", {})
        assert ns.get(required), f"{p.name}: nextstrain.{required} is unset"


def test_build_names_are_url_safe(configs):
    """build_name becomes a nextstrain.org path segment."""
    import re
    for p in configs:
        bn = yaml.safe_load(p.read_text())["nextstrain"]["build_name"]
        assert re.fullmatch(r"[A-Za-z0-9_-]+", bn), \
            f"{p.name}: build_name {bn!r} is not URL-safe"
