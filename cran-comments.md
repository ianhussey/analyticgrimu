## Submission

This is a new submission.

## Test environments

- local macOS, R release
- win-builder, R devel (via `devtools::check_win_devel()`)
- GitHub Actions: macOS, Windows, and Ubuntu (R devel, release, and oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note (new submission).

## Notes for the reviewer

* The package uses domain terms and author names that a spell-checker does
  not recognise but that are spelled correctly and intentional: GRIM, GRIM-U,
  Mann-Whitney, mid-p, StatXact, granularity, metascience, and the cited
  authors Heathers and Grimes. These are recorded in `inst/WORDLIST`, so
  `spelling::spell_check_package()` and the incoming spell-check pass with no
  spelling NOTE.
* The method reference (Heathers & Grimes, 2025) is a Medical Evidence
  Project report with no DOI; it is given in the Description as
  `Authors (year) <https://...>`.
* The package has no compiled code and writes no files outside `tempdir()`.
