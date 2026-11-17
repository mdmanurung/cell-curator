# annotation-consensus: determinism audit, 12 rows

Audited 2026-11-14 on HEAD `bc60f9e` (skill code unchanged since `15bed51`). Contract: `plans/v8-consensus-mode.md`, "Determinism contract". This audit replaces the earlier statement "passed 11 rows on the first snapshot", which `cc-v8-work/audit/fix-brief.md` contradicts (it lists defects on `snap-1239`).

Evidence type: hash comparison, no figure. `table fallback: the verdict is whether sha256 values are equal`.

## Verdict

Twelve of twelve rows pass for the demo fixture, with two limits stated below (rows 8 and 12).

| # | Requirement | Verdict | Evidence |
| --- | --- | --- | --- |
| 1 | Pinned environment; refuse drift unless `--allow-version-drift` | pass | tests `test_lock_file_matches_this_environment`, `test_version_drift_is_refused`, `test_allow_version_drift_is_recorded`, `test_verify_fails_when_the_lock_file_changes_after_the_run` |
| 2 | Every input hashed in `run_manifest.json` | pass | `test_manifest_holds_input_hashes`, `test_one_changed_byte_changes_the_manifest`, `test_verify_fails_when_an_input_changes` |
| 3 | No global RNG; independent of `PYTHONHASHSEED` | pass | source check `test_no_global_random_state`; new run: seeds 1 and 999 give the same `checksums.sha256` (sha256 `3ab725029cd9239a…`) and match all 18 files in `expected_hashes.json` |
| 4 | Identical at 1 and 8 threads | pass | `test_hashed_outputs_identical_at_threads_1_and_8`, `test_smoothed_votes_with_exact_duplicate_embeddings_ignore_row_order_and_threads` |
| 5 | Stable sorts, name alignment, near-tie rounding | pass | `test_row_order_of_inputs_does_not_matter`, kNN and tie tests |
| 6 | Fixed float formatting | pass | `test_every_float_in_every_hashed_table_has_at_most_the_configured_digits` |
| 7 | No time, host, user or path in hashed outputs | pass | `test_hashed_outputs_hold_no_host_time_or_temp_path`, `test_hashed_files_are_the_same_from_another_cwd_dir_name_user_and_config_location` |
| 8 | Login node and compute node give identical hashes | pass | job 26018992 on `res-hpc-exe024` against the login node `res-hpc-exe102`: same `checksums.sha256` hash, 18 of 18 files equal to `expected_hashes.json` |
| 9 | Snakemake and CLI give identical hashes | pass | `test_snakemake_and_cli_give_identical_hashed_outputs` (ran, 90 s, not skipped) |
| 10 | `verify` reruns in a scrubbed fresh process and fails on any difference | pass | `test_verify_reruns_in_a_subprocess_with_a_scrubbed_environment`, `test_verify_does_not_inherit_the_callers_environment`, `test_verify_fails_when_a_stored_output_file_is_edited` |
| 11 | Committed `expected_hashes.json` fails a changed rule | pass | `test_a_changed_rule_fails_the_expected_hashes` |
| 12 | No network at run time | pass, narrowed | run under `unshare -rn` (loopback only; a connection attempt gave "Network is unreachable"): exit 0, same hashes as the login run |

Fast suite on HEAD before this audit: 81 passed, 3 deselected, 336 s. Full suite earlier today: 83 passed, 1 skipped, 896 s.

## Container re-run (2026-11-16)

The image was rebuilt from `fd819fb` (`skills/annotation-consensus/container/`; SIF `2deafd07f26c…`). Rows 1, 8 and 12 were rerun inside it on the demo.

| Row | Inside the container | Evidence |
| --- | --- | --- |
| 1 | pass: the run starts with no `--allow-version-drift` and the manifest lists no differences | `container/_build/runs/*/run_manifest.json` |
| 8 | pass: login node, Xeon Gold 6234 (`res-hpc-path02`) and E5-4650 v3 (`res-hpc-exe032`) give 0 of 18 files differing from `expected_hashes.json` | `container/README.md` |
| 12 | pass: run under `unshare -rn`, 0 of 18 files differing | same |
| 10 | pass: `verify --demo` inside the container reports 18 stored files and the rerun byte-identical | same |

The earlier prototype differed on `environment.txt` (host 0.2.0 against wheel 0.3.0). That cause was stale install metadata in the host `.venv`. It was fixed in `fd819fb` (pin 0.3.0, `.venv` metadata refreshed), and all 18 files now match, `environment.txt` included.

### Rebuilt after the removal (2026-11-17)

The image was rebuilt from `a7083bb` (SIF `97febc139e1c…`, wheel `f6ddb0bd…`) because the package changed. Inside it the demo gives 0 of 18 differing files on the login node (AMD EPYC 7713P), under `unshare -rn` (row 12), and on `res-hpc-path01` and `res-hpc-path02` (Xeon Gold 6234). `verify --demo` passes (row 10) and the manifest holds the SIF hash. The run on `res-hpc-exe032` (E5-4650 v3, no AVX2) is not done: its job was queued 30+ minutes. Until it runs, row 8 inside the new image covers two CPU types, not three. Details: `skills/annotation-consensus/container/README.md`.

## Limits

- **Row 12 is narrower than the contract.** The contract lists CollecTRI, PROGENy and huCIRA resource files. The skill uses none of them: `annotation_consensus.py` never names them, and TF and pathway scores are not built. The run reads only the files in `input_checksums.sha256` (cells, votes, markers, tool predictions, vocabulary), and the manifest hashes each. The contract text still needs the user's wording decision.
- **Row 8 compares two nodes, not two CPU types.** I did not check the CPU model of `res-hpc-exe102` against `res-hpc-exe024`. A different BLAS kernel could still change a float. The planned container sets `OPENBLAS_CORETYPE` to cover this.
- **Rows 3, 8 and 12 ran the demo fixture only.** The AIDA and AIFI regressions (`ANNOTATION_CONSENSUS_FULL=1`) check scores against frozen numbers, not hashes across machines.
- **Older defects.** `cc-v8-work/audit/fix-brief.md` names defects in rows 1, 2, 4, 5, 6, 7 and 10 on `snap-1239`. Those rows pass on HEAD, but I did not trace each defect to its fix.

## Evidence files

`/exports/para-lipg-hpc/mdmanurung/cc-v8-work/audit/row8-2026-11-14/`: `compute.sbatch`, `logs/`, `out_compute/` (full output), and `checksums.sha256` for `out_login_seed1`, `out_login_seed999` and `out_nonet`.
