# Working with Forks, GitHub, and Codespaces (Student Guide)

> ✅ **TL;DR**
>
> 1. **Fork** the course repo.  
> 2. **Open Codespaces on your fork.**  
> 3. Create a **branch** (`feat/lab1-yourname`). Branch names always start with feat/<...> for featurees  
> 4. **Commit & push** to your fork.  
> 5. Open a **Pull Request** → base: `bozdogalex/bioinf-y4-lab:main`, head: `yourname/bioinf-y4-lab:lab1-yourname`.  
> 6. When the instructor updates `main`, **sync** via `git pull --ff-only upstream main` and **rebase your branch** if needed.

---

## 0) Before You Start

- Have a GitHub account and be logged in.
- Know your Git remotes vocabulary:  
  - **origin** = your fork (you can push to this)  
  - **upstream** = the instructor’s repo (you can only fetch/pull)

---

## 1) Fork the Repo

- Visit the course repo: `https://github.com/bozdogalex/bioinf-y4-lab`  
- Click **Fork** → **Create fork**.  
- Your fork will be: `https://github.com/<your-username>/bioinf-y4-lab`.

---

## 2) Open Codespaces **from Your Fork**

- Go to **your fork** (`<your-username>/bioinf-y4-lab`).  
- Click **Code → Codespaces → Create codespace on main**.

### How to verify you’re on your fork

- The **URL** should start with `github.com/<your-username>/...` (not `bozdogalex/...`), **or**
- In the terminal:
  ```bash
  git remote -v
  ```
  You should see:
  ```
  origin  https://github.com/<your-username>/bioinf-y4-lab.git
  ```
  (Later we’ll add `upstream` for the instructor’s repo.)

---

## 3) Set Up Remotes (one time)

Inside the Codespace terminal:

```bash
# Show current remotes
git remote -v

# If you don't see 'upstream', add the instructor's repo:
git remote add upstream https://github.com/bozdogalex/bioinf-y4-lab.git

# Verify
git remote -v
```

**Expected:**
```
origin   https://github.com/<your-username>/bioinf-y4-lab.git (fetch/push)
upstream https://github.com/bozdogalex/bioinf-y4-lab.git       (fetch)
```

> ⚠️ If `origin` shows `bozdogalex/...`, you opened Codespaces on the **wrong repo**. Close it and reopen Codespaces **from your fork**.  
> (Alternatively, you can `git remote rename origin upstream` then add your fork as `origin`, but it’s easier to start on the right repo.)

---

## 4) Create a Branch for Your Lab

```bash
git switch -c feat/lab1-<yourname>   # e.g., lab1-alex
```

Work, then commit & push:

```bash
git add -A
git commit -m "Lab 1: <short summary>"
git push -u origin HEAD
```

---

## 5) Open a Pull Request (PR) to the Instructor’s Repo

On **your fork** page → **Pull requests → New pull request** → click **Compare across forks** and set:

- **Base repository:** `bozdogalex/bioinf-y4-lab`  
- **Base branch:** `main` *(unless told to use a lab-specific base branch)*  
- **Head repository:** `<your-username>/bioinf-y4-lab`  
- **Head branch:** `lab1-<yourname>`

Then **Create pull request**.

**PR title suggestion:**
```
lab1(<yourname>): short summary
```

---

## 6) Keep Your Fork Up to Date (Sync with Instructor)

Run this regularly, especially after the instructor merges fixes:

```bash
# Update your local 'main' from the instructor's 'main'
git fetch upstream
git switch main
git pull --ff-only upstream main

# Push the updated main back to your fork
git push origin main
```

### Update your lab branch onto the new main (when PR shows “out of date”)

```bash
git switch lab1-<yourname>
git rebase main
# Resolve any conflicts, then:
git push --force-with-lease
```

✅ `--ff-only` on `main` keeps history clean.  
✅ Use **rebase** for your feature/lab branches when needed.

### Web UI alternative

On your fork’s GitHub page → **Sync fork / Fetch upstream → Update branch**.  
Then in Codespaces:

```bash
git pull
```

---

## 7) When the Instructor Updates Your Open PR Branch

If the PR allows **“edits by maintainers,”** the instructor may push changes to your PR branch.

- Codespaces may prompt you to **Pull**.  
- Or pull manually:
  ```bash
  git pull
  ```

No need to reopen Codespaces.

---

## 8) Environment Changes (requirements/devcontainer)

If `requirements.txt`, `.devcontainer/`, or `Dockerfile` changed:

1. Pull the latest code (Section 6).  
2. Rebuild the container (no need to recreate Codespaces):  
   **Ctrl+Shift+P** → “**Codespaces: Rebuild Container**”.

---

## 9) Common Errors & Fixes

### “origin does not seem to be a repository”

You don’t have `origin` set, or it’s wrong.

```bash
git remote -v
# If missing:
git remote add origin https://github.com/<your-username>/bioinf-y4-lab.git
```

### `git pull` says “Already up to date” but instructor merged changes

You’re pulling from **origin** (your fork), not from **upstream**.

```bash
git fetch upstream
git switch main
git pull --ff-only upstream main
git push origin main
```

### Pushed to wrong place / opened Codespaces on instructor repo

Check:

```bash
git remote -v
```

If `origin` shows `bozdogalex/...`, you’re on the wrong repo.  
Best fix: open a new Codespace **from your fork**.

### “Permission denied (publickey)” with SSH

Use HTTPS URLs in remotes, or set up SSH keys first.

### Merge conflicts during rebase

Git will mark conflicts. Fix files, then:

```bash
git add <fixed-file>
git rebase --continue
```

If you need to abort:

```bash
git rebase --abort
```

---

## 10) Minimal Daily Commands (Cheat Sheet)

```bash
# Sync my main with instructor's main
git fetch upstream
git switch main
git pull --ff-only upstream main
git push origin main

# Update my lab branch to the latest main
git switch lab1-<yourname>
git rebase main
git push --force-with-lease

# Create PR from my branch to instructor/main (via GitHub UI)
```

---

If anything here is unclear, ask on the course discussion board with your PR link and the output of `git remote -v`.
