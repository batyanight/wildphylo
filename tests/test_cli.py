"""The ``wildphylo`` command dispatches to registered stage steps."""

import argparse

import pytest

from wildphylo import cli


def _add_arguments_echo(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--value", default="x")


def _run_echo(args: argparse.Namespace) -> int:
    print(f"echo {args.value}")
    return 0


def _run_fail(_args: argparse.Namespace) -> int:
    return 2


@pytest.fixture
def registry(monkeypatch):
    """A clean registry with one passing and one failing step under 'gates'."""
    fresh = {stage: {} for stage in cli.STAGES}
    monkeypatch.setattr(cli, "_REGISTRY", fresh)
    cli.register(cli.Step("gates", "echo", "echo the value", _add_arguments_echo, _run_echo))
    cli.register(cli.Step("gates", "fail", "always fail", _add_arguments_echo, _run_fail))
    return fresh


def test_every_stage_in_the_plan_has_a_subparser():
    parser = cli.build_parser()
    stage_action = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))  # noqa: SLF001
    assert set(stage_action.choices) == set(cli.STAGES)
    assert cli.STAGES[:8] == (
        "fetch",
        "curate",
        "align",
        "tree",
        "gates",
        "beast",
        "auspice",
        "publish",
    )


@pytest.mark.usefixtures("registry")
def test_a_step_runs_with_its_own_arguments(capsys):
    assert cli.main(["gates", "echo", "--value", "hello"]) == 0
    assert capsys.readouterr().out.strip() == "echo hello"


@pytest.mark.usefixtures("registry")
def test_a_failing_gate_returns_its_exit_code():
    assert cli.main(["gates", "fail"]) == 2


@pytest.mark.usefixtures("registry")
def test_an_unknown_step_is_a_usage_error():
    with pytest.raises(SystemExit) as exc:
        cli.main(["gates", "nope"])
    assert exc.value.code == 2


@pytest.mark.usefixtures("registry")
def test_duplicate_registration_is_an_error():
    with pytest.raises(ValueError, match="registered twice"):
        cli.register(cli.Step("gates", "echo", "again", _add_arguments_echo, _run_echo))


@pytest.mark.usefixtures("registry")
def test_the_decorator_needs_an_add_arguments_hook():
    with pytest.raises(AttributeError, match="add_arguments_orphan"):

        @cli.step("gates", "orphan", "no hook")
        def _orphan(_args: argparse.Namespace) -> int:
            return 0
