"""
The shared ``--build-dir`` mechanics: flags, precedence, validation, usage errors.

The per-step tables (which flag comes from which layout method) are tested with
each step. Here only the helper itself is, on a toy step, so a failure points at
``wildphylo.stepargs`` and not at one stage.
"""

import argparse
import logging
from pathlib import Path

import pytest

from wildphylo.layout import BuildLayout, dataset_path
from wildphylo.stepargs import (
    Build,
    UsageError,
    add_build_dir_arguments,
    fail,
    resolve,
)

CONFIG = Path("tests/fixtures/configs/cdv.yaml")


def parse(*argv: str, **kwargs: bool) -> argparse.Namespace:
    """A toy step: ``--tree``, ``--out`` and ``--n`` plus the shared flags."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--n", type=int)
    add_build_dir_arguments(parser, **kwargs)
    return parser.parse_args(argv)


def tree_table() -> dict:
    return {
        "tree": lambda b: b.layout.ml_tree(b.need("locus")),
        "out": lambda b: b.layout.temporal_report(b.need("locus")),
        "n": lambda b: b.setting("temporal_signal", "permutations"),
    }


# --- add_build_dir_arguments ------------------------------------------------


def test_default_flags_are_build_dir_config_and_locus():
    args = parse("--build-dir", "b", "--config", "c.yaml", "--locus", "H")
    assert (args.build_dir, args.config, args.locus) == (Path("b"), Path("c.yaml"), "H")
    assert not hasattr(args, "seed")
    assert not hasattr(args, "rep")


def test_flags_default_to_none():
    args = parse()
    assert (args.build_dir, args.config, args.locus) == (None, None, None)


def test_seed_is_an_int_and_rep_a_string():
    args = parse("--seed", "7", "--rep", "03", seed=True, rep=True)
    assert args.seed == 7
    assert args.rep == "03"


def test_locus_can_be_left_out():
    with pytest.raises(SystemExit):
        parse("--locus", "H", locus=False)


def test_an_existing_config_flag_is_kept():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True, help="mine")
    add_build_dir_arguments(parser)  # must not clash
    assert parser.parse_args(["--config", "c"]).config == Path("c")
    assert "mine" in parser.format_help()


# --- precedence -------------------------------------------------------------


def test_without_build_dir_nothing_is_derived():
    args = parse("--tree", "t.nwk")
    assert resolve(args, derived=tree_table()) is None
    assert args.out is None
    assert args.n is None


def test_derived_values_fill_what_is_missing():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--locus", "H")
    build = resolve(args, derived=tree_table())
    assert isinstance(build, Build)
    assert args.tree == build.layout.ml_tree("H")
    assert args.out == build.layout.temporal_report("H")
    assert args.n == build.setting("temporal_signal", "permutations")
    assert build.layout == BuildLayout.at("b", build.cfg)


def test_explicit_flags_beat_derived_values_and_are_not_computed():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--tree", "mine.nwk", "--n", "3")

    def boom(_: Build) -> str:
        raise AssertionError("an explicit flag must not call its default")

    resolve(args, derived={**tree_table(), "tree": boom, "n": boom, "out": lambda _: Path("o")})
    assert args.tree == Path("mine.nwk")
    assert args.n == 3  # also without --locus: nothing needed it


def test_a_falsy_explicit_value_still_wins():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--n", "0")
    resolve(args, derived={"n": lambda _: 99}, fallbacks={"n": 5})
    assert args.n == 0


def test_fallbacks_rank_below_derived_values_and_above_nothing():
    args = parse("--build-dir", "b", "--config", str(CONFIG))
    resolve(args, derived={"n": lambda _: 4}, fallbacks={"n": 1, "out": Path("legacy.json")})
    assert args.n == 4
    assert args.out == Path("legacy.json")


def test_a_derived_none_means_no_default_so_the_fallback_applies():
    args = parse("--build-dir", "b", "--config", str(CONFIG))
    resolve(args, derived={"n": lambda b: b.setting("no", "such", "key")}, fallbacks={"n": 1000})
    assert args.n == 1000


def test_without_build_dir_fallbacks_still_apply():
    args = parse()
    resolve(args, derived=tree_table(), fallbacks={"n": 1000})
    assert args.n == 1000


def test_a_callable_fallback_sees_the_values_resolved_before_it():
    args = parse("--out", "x/result.json")
    resolve(
        args,
        fallbacks={"tree": lambda a: a.out.with_suffix(".nwk"), "n": lambda a: len(str(a.tree))},
    )
    assert args.tree == Path("x/result.nwk")
    assert args.n == len("x/result.nwk")


def test_a_derived_value_may_use_one_resolved_before_it():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--locus", "H")
    resolve(
        args,
        derived={
            "tree": lambda b: b.layout.ml_tree(b.need("locus")),
            "out": lambda _: args.tree.with_suffix(".json"),
        },
    )
    assert args.out == args.tree.with_suffix(".json")


# --- validation -------------------------------------------------------------


def test_build_dir_needs_config():
    args = parse("--build-dir", "b", "--locus", "H")
    with pytest.raises(UsageError, match="--build-dir needs --config"):
        resolve(args, derived=tree_table())


def test_a_missing_config_file_is_a_usage_error():
    args = parse("--build-dir", "b", "--config", "nope.yaml")
    with pytest.raises(UsageError, match="config not found"):
        resolve(args, derived=tree_table())


def test_a_preloaded_config_replaces_the_config_flag():
    args = parse("--build-dir", "b", "--locus", "H")
    build = resolve(args, derived=tree_table(), cfg={"id": "toy"})
    assert build is not None
    assert build.layout.pathogen == "toy"


def test_the_config_is_validated_not_taken_as_written(tmp_path):
    config = tmp_path / "bare.yaml"
    config.write_text("id: bare\n")  # not a complete pathogen config
    args = parse("--build-dir", "b", "--config", str(config), "--locus", "H")
    with pytest.raises(UsageError, match="not usable"):
        resolve(args, derived={"tree": lambda b: b.layout.ml_tree("H")})


def test_a_missing_wildcard_flag_is_named_with_everything_that_needed_it():
    args = parse("--build-dir", "b", "--config", str(CONFIG))
    with pytest.raises(UsageError) as err:
        resolve(args, derived=tree_table())
    assert str(err.value) == "--build-dir needs --locus to derive --tree, --out (or give them)"


def test_a_wildcard_is_not_needed_when_every_path_is_given():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--tree", "t", "--out", "o")
    resolve(args, derived={**tree_table(), "n": lambda _: 1})
    assert args.n == 1


def test_seed_and_rep_are_wildcards_too():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--locus", "H", seed=True, rep=True)
    derived = {"out": lambda b: b.layout.beast_xml(b.need("locus"), b.need("seed"))}
    with pytest.raises(UsageError, match="--seed to derive --out"):
        resolve(args, derived=derived)
    args.seed = 0  # a zero seed is a seed
    resolve(args, derived=derived)
    assert args.out == Path("b/beast/H/seed0/cdv.xml")


def test_an_error_inside_a_default_names_the_flag_it_was_for():
    def broken(_: Build) -> Path:
        msg = "the config has no such thing"
        raise ValueError(msg)

    args = parse("--build-dir", "b", "--config", str(CONFIG))
    with pytest.raises(UsageError, match="cannot default --tree from --build-dir: the config has"):
        resolve(args, derived={"tree": broken})


def test_required_flags_are_reported_together():
    args = parse()
    with pytest.raises(UsageError) as err:
        resolve(args, required=("tree", "out"))
    assert str(err.value) == "--tree, --out required without --build-dir"


def test_one_missing_required_flag_reads_as_a_sentence():
    args = parse("--out", "o")
    with pytest.raises(UsageError) as err:
        resolve(args, required=("tree", "out"))
    assert str(err.value) == "--tree is required without --build-dir"


def test_required_is_met_by_a_derived_value():
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--locus", "H")
    resolve(args, derived=tree_table(), required=("tree", "out"))


def test_locus_only_applies_with_build_dir_when_the_step_says_so():
    args = parse("--locus", "H")
    with pytest.raises(UsageError, match="--locus only applies together with --build-dir"):
        resolve(args, build_dir_only=("locus",))
    resolve(args)  # other steps use --locus on its own


# --- outputs ----------------------------------------------------------------


def test_derived_outputs_get_their_parent_directories(tmp_path):
    args = parse("--build-dir", str(tmp_path / "b"), "--config", str(CONFIG), "--locus", "H")
    resolve(args, derived=tree_table(), outputs=("out",))
    assert args.out.parent.is_dir()
    assert not args.tree.parent.exists()  # an input: not ours to create


def test_explicit_outputs_are_left_alone(tmp_path):
    out = tmp_path / "elsewhere" / "o.json"
    build_dir = str(tmp_path / "b")
    args = parse(
        "--build-dir", build_dir, "--config", str(CONFIG), "--locus", "H", "--out", str(out)
    )
    resolve(args, derived=tree_table(), outputs=("out",))
    assert not out.parent.exists()


def test_nothing_is_created_when_the_command_line_is_refused(tmp_path):
    args = parse("--build-dir", str(tmp_path / "b"), "--config", str(CONFIG), "--locus", "H")
    with pytest.raises(UsageError):
        resolve(
            args,
            derived={**tree_table(), "n": lambda _: None},
            outputs=("out",),
            required=("n",),
        )
    assert not (tmp_path / "b").exists()


# --- logging and reporting --------------------------------------------------


def test_each_resolved_value_is_logged_at_debug_level(caplog):
    args = parse("--build-dir", "b", "--config", str(CONFIG), "--locus", "H")
    with caplog.at_level(logging.DEBUG, logger="wildphylo.stepargs"):
        resolve(args, derived=tree_table(), fallbacks={"n": 1})
    text = caplog.text
    assert "--tree = b/tree/H_ml.treefile (from --build-dir)" in text
    assert "--out = " in text
    assert all(r.levelno == logging.DEBUG for r in caplog.records)


def test_fail_logs_the_message_and_returns_the_code(caplog):
    with caplog.at_level(logging.ERROR, logger="wildphylo.stepargs"):
        assert fail(UsageError("--tree is required")) == 2
    assert "--tree is required" in caplog.text
    assert [r.levelno for r in caplog.records] == [logging.ERROR]


# --- Build ------------------------------------------------------------------


def make_build(cfg: dict, **flags: object) -> Build:
    return Build(BuildLayout.at("b", cfg), cfg, argparse.Namespace(**flags))


def test_need_returns_the_flag_as_a_string_and_refuses_none_and_empty():
    build = make_build({"id": "x"}, locus="H", seed=0, rep=None, empty="")
    assert build.need("locus") == "H"
    assert build.need("seed") == "0"
    for flag in ("rep", "empty", "absent"):
        with pytest.raises(UsageError, match=f"--{flag} is needed"):
            build.need(flag)


def test_setting_walks_the_config_and_gives_none_off_the_path():
    build = make_build({"id": "x", "a": {"b": 3, "c": "leaf"}})
    assert build.setting("a", "b") == 3
    assert build.setting("a", "zzz") is None
    assert build.setting("a", "c", "deeper") is None
    assert build.setting() == build.cfg


def test_loci_are_the_analysed_ones_or_the_single_locus():
    assert make_build({"id": "x", "locus": {"name": "H"}}).loci() == ["H"]
    cfg = {"id": "x", "loci": [{"name": "L", "analyse": True}, {"name": "M"}]}
    assert make_build(cfg).loci() == ["L"]


def test_seeds_are_strings_and_an_empty_list_is_an_error():
    assert make_build({"id": "x", "beast": {"seeds": [1, 2]}}).seeds() == ["1", "2"]
    with pytest.raises(UsageError, match=r"beast.seeds is empty"):
        make_build({"id": "x"}).seeds()


# --- dataset_path -----------------------------------------------------------


def test_a_single_locus_keeps_the_bare_dataset_name():
    cfg = {"id": "x", "nextstrain": {"build_name": "x-phylo"}, "locus": {"name": "H"}}
    layout = BuildLayout.at("b", cfg)
    assert dataset_path(layout, cfg, "H") == Path("b/auspice/x-phylo.json")
    assert dataset_path(layout, cfg, None) == Path("b/auspice/x-phylo.json")


def test_several_analysed_loci_get_one_dataset_each():
    cfg = {
        "id": "x",
        "nextstrain": {"build_name": "x-phylo"},
        "loci": [{"name": "L", "analyse": True}, {"name": "M", "analyse": True}, {"name": "S"}],
    }
    layout = BuildLayout.at("b", cfg)
    assert dataset_path(layout, cfg, "L") == Path("b/auspice/x-phylo_L.json")
    with pytest.raises(ValueError, match="locus is needed"):
        dataset_path(layout, cfg, None)


def test_a_config_without_a_build_name_cannot_name_the_dataset():
    cfg = {"id": "x", "locus": {"name": "H"}}
    with pytest.raises(ValueError, match=r"nextstrain.build_name"):
        dataset_path(BuildLayout.at("b", cfg), cfg, "H")
