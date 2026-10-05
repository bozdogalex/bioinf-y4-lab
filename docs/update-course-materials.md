# Updating an existing private coursework repository

Adding the instructor as a collaborator grants access; it does not synchronize repositories. Keep your existing private repository and invitation. Do not create another repository or make it public to receive updates.

The private repositories were created from snapshots, so their histories may be unrelated to the public course history. We import an explicit list of teaching files on a new private branch. We do not merge public history or restore whole lab directories.

## Before starting

Save and commit your current work **on its current branch**, reviewing `git status` first. Keep working datasets in the ignored data directory. If you edited an original teaching script, preserve that work in your personal submission directory before updating: listed teaching files will be replaced, while personal submissions and roster entries are excluded.

In your private Codespace, confirm that `git remote -v` shows your private repository as `origin`. The instructor supplies a tested course commit SHA. After the update has been merged into the public repository, `main` can also be used. Until then, use the instructor-supplied SHA; the old public `main` does not contain this update.

## Import the Lab 1 and 2 update

Set `COURSE_REF` to the instructor's commit SHA (or `main` after its release). Then paste the following Bash block into your Codespaces terminal. It stops on an error and refuses an uncommitted working tree. If the update branch already exists, resume it rather than deleting it.

```bash
read -r -p "Course commit SHA (or main after release): " COURSE_REF
export COURSE_REF
(
  set -eu
  cd "$(git rev-parse --show-toplevel)"
  test -n "$COURSE_REF"
  if test -n "$(git status --porcelain)"; then
    echo "STOP: commit or otherwise preserve your changes on the current branch first."
    exit 1
  fi
  git switch main
  git pull --ff-only origin main
  git switch -c update/lab12-2026-10-05
  git fetch --no-tags --depth=1 https://github.com/bozdogalex/bioinf-y4-lab.git "$COURSE_REF"
  COURSE_SHA=$(git rev-parse FETCH_HEAD)
  UPDATE_DIR=$(mktemp -d)
  UPDATE_LIST="$UPDATE_DIR/lab12.paths"
  git show "$COURSE_SHA:docs/updates/lab12-2026-10-05.paths" > "$UPDATE_LIST"
  python - "$UPDATE_LIST" <<'PY'
from pathlib import Path, PurePosixPath
import sys
paths = Path(sys.argv[1]).read_text().splitlines()
if not paths or len(paths) != len(set(paths)):
    raise SystemExit("Invalid or duplicate update paths")
for item in paths:
    parts = PurePosixPath(item).parts
    if (not parts or PurePosixPath(item).is_absolute() or ".." in parts
            or "\\" in item or ":" in item
            or any(p in {".git", "submissions", "roster"} for p in parts)
            or item.startswith("data/work/")):
        raise SystemExit(f"Refusing protected/unsafe path: {item}")
PY
  git --literal-pathspecs restore --source="$COURSE_SHA" --worktree --pathspec-from-file="$UPDATE_LIST"
  git --literal-pathspecs add --pathspec-from-file="$UPDATE_LIST"
  git diff --cached --stat
  git diff --cached --name-only
  echo "Imported course commit: $COURSE_SHA"
  echo "Review the staged changes, run the checks, then commit and push this private branch."
)
```

The file list names individual scripts, corrected datasets, documentation and checks. It excludes all `submissions/`, `roster/` and `data/work/` contents. A prior version of any replaced tracked file remains in your private Git history. The command does not rewrite history or push anywhere.

## Review and finish

```bash
git diff --cached
python labs/00_smoke/smoke.py
python -m unittest discover -s tests -p 'test_lab12*.py' -v
python "labs/01_intro&databases/demo01_entrez_brca1.py"
python labs/02_alignment/demo01_pairwise.py --fasta data/sample/toy_alignment.fasta --k 10
git commit -m "Update Lab 1 and 2 teaching materials"
git push -u origin HEAD
```

Create an update PR **inside your private repository**, targeting its `main`, review it and merge the teaching update. This is a materials update, not an assessed lab submission. Do not include or merge pending solutions merely to receive the update.

After merging that update PR:

```bash
git switch main
git pull --ff-only origin main
```

For an existing exercise branch, switch back to its actual name and run `git merge main`. This merges your own private main, not the public course history. It brings the teaching update into the exercise branch without a force push. If there is a conflict, stop and ask for help; do not discard your solution with a hard reset. The update does not change the container dependencies, so rebuilding Codespaces is not required for this release.

## If you have not finished or committed your work

Do not use `git reset --hard`, replace the entire lab folder, or extract a ZIP over your coursework. Preserve your work first. `data/work/` is ignored by Git and is not included in a normal commit; the update deliberately leaves it untouched. Keep any important independent backup required by your normal workflow.

For later releases, use the new manifest named by the instructor. This manifest applies only to the Lab 1 and 2 update of 5 October 2026.
