# Private student repositories — 2026–2027

The public course repository contains teaching materials. Assessed work stays in your own **private repository**, shared with the instructor **bozdogalex**. Do not submit solutions or your roster entry to the public course repository.

## 1. Create your private repository

1. On the course repository, choose **Code → Download ZIP** and extract it. Use the cleaned 2026–2027 version.
2. Create a new GitHub repository named `bioinf-y4-2026-2027-<handle>` with **Private** visibility. Do not initialize it with a README, license, or .gitignore.
3. Open a terminal in the extracted folder containing the course README. Initialize and upload this snapshot (replace `<handle>`):

```bash
git init -b main
git add .
git commit -m "Initialize private BIOINF-Y4 coursework"
git remote add origin https://github.com/<handle>/bioinf-y4-2026-2027-<handle>.git
git push -u origin main
```

Create a new repository, not a public fork. The ZIP contains the current files without the course Git history. Preserve the included configuration files, license, and attribution. If a PDF is only a Git LFS pointer in your download, open/download the actual file from the public course repository.

4. In your private repository, open **Settings → Collaborators**, invite **bozdogalex**, and ensure the invitation is accepted before assessment.
5. Submit the private repository URL through the university LMS. Keep the repository private throughout the course.

## 2. Open your working environment

Open **Code → Codespaces → Create codespace on main** from your private repository, or clone it locally and use Docker. See [onboarding](onboarding.md). Check `git remote -v`: `origin` must point to your private repository.

The copied CI workflow can check your work if GitHub Actions is enabled and your account has available usage. Otherwise run the same checks locally and include the result in your PR. The image-publishing workflow runs only in the instructor's repository; students use the existing course image.

## 3. One branch and PR per lab

Start each lab from your private `main`:

```bash
git switch main
git pull --ff-only origin main
git switch -c feat/lab02-<handle>
```

Copy the exercise skeleton into `labs/NN_topic/submissions/<handle>/` and complete it there. Keep the teaching skeleton unchanged. Follow the lab's deliverable requirements; use `data/work/<handle>/` for local working datasets and do not commit large or sensitive data.

For the first lab, add your own row to `labs/01_intro&databases/roster/handles.csv` **only in your private repository**.

```bash
git add labs
git commit -m "Lab 02: submission"
git push -u origin HEAD
```

Open a PR **within your private repository**:

- Base repository: your private repository; base branch: `main`.
- Head repository: the same private repository; head branch: `feat/lab02-<handle>`.
- Title: `Lab 02 — <handle>`.
- Fill in the PR checklist and describe how to reproduce your results.

Do not select `bozdogalex/bioinf-y4-lab` as the base repository for coursework.

## 4. Submit for assessment

By the deadline announced in the LMS, submit:

- the private PR URL;
- the full submitted commit SHA (`git rev-parse HEAD`);
- any additional items explicitly requested by the instructor.

Keep the PR open for review. The instructor can inspect code, results, CI, and the recorded deadline commit, and request corrections. Push corrections to the same branch; do not rewrite submitted history. Grades are recorded in the LMS. A green CI result checks the environment/syntax, not the scientific correctness of your answer.

Merge into your private `main` after assessment or when the instructor allows it. If the next lab depends on work still awaiting review, branch from that lab's branch and explain the dependency in the next PR. Do not wait for grading to continue your work.

## 5. Receive teaching updates

For the 5 October 2026 Lab 1/2 corrections, follow the [tested selective update procedure](update-course-materials.md). It imports a named list of teaching files on a new branch while excluding submissions, roster entries and working data. Existing collaborator invitations remain valid. Use the instructor-supplied commit SHA until the update is released on public main.

Check the public course repository for announcements and updated files. Download the updated materials and copy only the files the instructor identifies into your private repository, reviewing changes before committing. Preserve your submissions and personal roster row. Do not merge the public repository's historical branches into your private repository.

## Instructor assessment

Use the LMS repository/PR links to access each private submission. Review the recorded commit for deadline assessment, leave feedback on the private PR, and record the grade in the LMS. Never merge student solutions into the public teaching repository. Public PRs are for teaching-material fixes and improvements only.
