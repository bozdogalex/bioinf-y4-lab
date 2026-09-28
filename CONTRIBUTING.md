# Contributing Guidelines — BIOINF-Y4 Lab

## Coursework (2026–2027)

Complete assessed work in your own **private repository**, shared with the instructor `bozdogalex`. Open each lab PR **inside that private repository**, targeting its `main` branch. Submit the PR URL and deadline commit SHA through the university LMS.

Follow the [private repository and submission guide](docs/git-workflow.md) for setup, instructor access, review, and receiving teaching updates. Do not send solutions, reports, or roster entries to the public course repository.

- Put completed exercise copies and results in `labs/NN_topic/submissions/<handle>/`.
- Preserve the original teaching skeletons.
- Follow each lab's deliverable requirements; reports are at most two pages unless specified otherwise.
- Keep large working datasets in `data/work/<handle>/`; do not upload sensitive data.
- Include execution instructions, results, and any AI-assistance attribution.
- Run syntax, smoke, and MLflow checks where available; CI is not a correctness grade.
- Keep assessment PRs open until reviewed or instructed otherwise. Grades are recorded in the LMS.

## Pair work

Labs may be completed individually or in pairs as allowed by the instructor. For pairs, identify both contributors and rotate driver/navigator roles. Confirm the submission arrangement with the instructor; do not expose a partner's work publicly.

## Public teaching-material contributions

Public PRs are welcome for corrections and improvements to teaching material. Fork the public repository, use a descriptive branch, and explain the change. Do not include assessed solutions, generated student results, personal roster data, or grades. Maintainers review teaching changes before merging.

## Style and data policy

Use readable Python and clear commit messages. Clear notebook outputs before committing unless required for assessment. Follow [data policy](docs/GDPR_and_DataPolicy.md) and [repository policies](docs/policies.md).
