#!/usr/bin/env bash
#
# Create the GitHub repository and push this scaffold.
#
#   bash INIT_REPO.sh --check              inspect only, change nothing
#   bash INIT_REPO.sh                      create as "wildphylo"
#   bash INIT_REPO.sh myreponame           create under a different name
#   bash INIT_REPO.sh myreponame --private
#
# Safe to re-run. Each step checks whether it has already been done, so an
# interrupted run can be resumed rather than restarted.

set -euo pipefail

REPO=""
VISIBILITY="--public"
CHECK_ONLY=0

for arg in "$@"; do
  case "$arg" in
    --check)    CHECK_ONLY=1 ;;
    --private)  VISIBILITY="--private" ;;
    --public)   VISIBILITY="--public" ;;
    -*)         echo "unknown option: $arg" >&2; exit 2 ;;
    *)          REPO="$arg" ;;
  esac
done
REPO="${REPO:-wildphylo}"

say()  { printf '\n\033[1m%s\033[0m\n' "$*"; }
ok()   { printf '  \033[32mok\033[0m    %s\n' "$*"; }
warn() { printf '  \033[33mwarn\033[0m  %s\n' "$*"; }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$*"; }

# ---------------------------------------------------------------------------
# Preflight. Everything that could fail mid-run is checked up front, because a
# half-created repo is more annoying to clean up than a refused one.
# ---------------------------------------------------------------------------
say "Preflight"
FAIL=0
OWNER=""

if [[ ! -f workflow/Snakefile || ! -d config/pathogen ]]; then
  bad "run this from the repository root (workflow/ and config/ must be here)"
  FAIL=1
else
  ok "in the repository root"
fi

for tool in git gh; do
  if command -v "$tool" >/dev/null 2>&1; then
    ok "$tool found"
  else
    bad "$tool not found"
    [[ $tool == gh ]] && echo "        install: https://cli.github.com/  (brew install gh)"
    FAIL=1
  fi
done

if command -v gh >/dev/null 2>&1; then
  if gh auth status >/dev/null 2>&1; then
    OWNER="$(gh api user --jq .login 2>/dev/null || echo '')"
    ok "gh authenticated as ${OWNER:-unknown}"
  else
    bad "gh not authenticated — run:  gh auth login"
    FAIL=1
  fi
fi

if git rev-parse --git-dir >/dev/null 2>&1; then
  warn "already a git repository; existing history will be kept"
  GIT_EXISTS=1
else
  ok "not yet a git repository"
  GIT_EXISTS=0
fi

# An unset identity aborts the commit AFTER git init has run, which leaves a
# half-initialised repo. Check it here instead.
GIT_NAME="$(git config --get user.name  || true)"
GIT_MAIL="$(git config --get user.email || true)"
if [[ -n "$GIT_NAME" && -n "$GIT_MAIL" ]]; then
  ok "git identity: $GIT_NAME <$GIT_MAIL>"
else
  bad "git identity not set — run:"
  echo "        git config --global user.name  \"Your Name\""
  echo "        git config --global user.email \"you@example.com\""
  FAIL=1
fi

REPO_EXISTS=0
if [[ -n "$OWNER" ]]; then
  if gh repo view "$OWNER/$REPO" >/dev/null 2>&1; then
    warn "$OWNER/$REPO already exists; will push to it rather than recreate"
    REPO_EXISTS=1
  else
    ok "$OWNER/$REPO is available"
  fi
fi

# Run the tests before publishing. A scaffold whose own suite fails is a bad
# first commit to have permanently in the history.
if command -v pytest >/dev/null 2>&1; then
  if pytest -q >/dev/null 2>&1; then
    ok "test suite passes"
  else
    bad "test suite fails — fix before publishing (run: pytest)"
    FAIL=1
  fi
else
  warn "pytest not found; skipping the test check (pip install -e '.[dev]')"
fi

if [[ $FAIL -eq 1 ]]; then
  say "Preflight failed. Nothing was changed."
  exit 1
fi

if [[ $CHECK_ONLY -eq 1 ]]; then
  say "Check only — nothing was changed."
  echo "  Would create: ${OWNER:-<you>}/$REPO  ($VISIBILITY)"
  exit 0
fi

# ---------------------------------------------------------------------------
say "Committing"
# ---------------------------------------------------------------------------
if [[ $GIT_EXISTS -eq 0 ]]; then
  git init -b main >/dev/null
  ok "initialised"
fi

# git init -b needs git >= 2.28; rename if the flag was ignored.
git symbolic-ref -q HEAD refs/heads/main >/dev/null 2>&1 || git branch -M main

git add -A
if git diff --staged --quiet; then
  ok "nothing to commit (already committed)"
else
  git commit -q -m "Initial commit: config-driven wildlife phylodynamics pipeline

Extracted and generalised from cdv-phylodynamics.

- Config-driven: every pathogen constant lives in config/pathogen/<id>.yaml
- Runtime gates for manual review, temporal signal, segment congruence,
  MCMC convergence, and build-to-build comparison
- 144 tests, including regression tests for three bugs found in the
  original implementation
- Two worked configs: CDV (unsegmented) and BTV (10 segments)"
  ok "committed"
fi

# ---------------------------------------------------------------------------
say "Publishing"
# ---------------------------------------------------------------------------
if [[ $REPO_EXISTS -eq 0 ]]; then
  gh repo create "$REPO" $VISIBILITY --source=. --remote=origin \
    --description "Config-driven phylodynamics pipeline for wildlife disease" \
    --push
  ok "created and pushed"
else
  git remote get-url origin >/dev/null 2>&1 || \
    git remote add origin "https://github.com/$OWNER/$REPO.git"
  git push -u origin main
  ok "pushed to the existing repository"
fi

gh repo edit "$REPO" \
  --add-topic phylodynamics --add-topic wildlife-disease \
  --add-topic nextstrain --add-topic beast \
  --add-topic molecular-epidemiology --add-topic snakemake \
  --add-topic reproducible-research >/dev/null
ok "topics set"

# ---------------------------------------------------------------------------
say "Next steps (manual)"
# ---------------------------------------------------------------------------
cat <<EOT

  1. SECRETS — needed by the scheduled rebuild workflow:

       gh secret set NCBI_EMAIL    --body "you@example.com"
       gh secret set NCBI_API_KEY  --body "YOUR_KEY"

     An API key is free and raises the Entrez rate limit from 3/s to 10/s:
     https://www.ncbi.nlm.nih.gov/account/
     Never put these in a config file.

  2. BRANCH PROTECTION — deliberately NOT set by this script.

     Requiring a status check that never runs blocks every merge permanently,
     and the exact check name depends on the CI matrix. Let CI run once, then
     copy the real job names from the Actions tab into the "contexts" list:

       gh api -X PUT repos/$OWNER/$REPO/branches/main/protection --input - <<'JSON'
       {
         "required_status_checks": {
           "strict": true,
           "contexts": ["tests (py3.13)", "validate shipped configs"]
         },
         "enforce_admins": false,
         "required_pull_request_reviews": null,
         "restrictions": null
       }
       JSON

     Skip this entirely if you are working solo; it mostly gets in the way.

  3. ZENODO — link the repo at
     https://zenodo.org/account/settings/github/
     then cut a release to mint a DOI:

       gh release create v0.1.0 --generate-notes

  4. LEAVE cdv-phylodynamics ALONE for now. Its Nextstrain Community URL is
     already published and should keep working. Strip the pipeline code out of
     it only once this repo is installable.

  Repository: https://github.com/$OWNER/$REPO

EOT
