# Cover letter — required edits

Two files are in play:

- `_tools/cover_letter_SR_v2.md` — the markdown source, **already corrected** by this
  pass.
- `Cover letter_Scientific Reports_修订版.docx` — the version dated 29 September 2026,
  which is what would actually be sent. **Still contains the three problems below.**

Apply the same three edits to the `.docx` before submission. Replacement text is given
verbatim.

---

## ⛔ 1. The letter asserts a public repository and a DOI that do not exist yet

**Current text (third paragraph of the letter):**

> **Data and code availability.** All datasets are publicly available under the GEO
> accessions and the eQTL and GWAS resources listed in the Data availability statement.
> The analysis code, including pre-specified parameters and random seeds, **has been
> deposited in a public repository, and the statement provides the permanent DOI.** *(If
> deposition is not yet complete, please substitute: The analysis code, including
> pre-specified parameters and random seeds, will be deposited in a public repository, and
> the Data availability statement provides the permanent link.)*

**The problem.** The repository is prepared locally but has never been pushed, is not
public, and has no DOI. As written, the letter tells the editorial office that both exist.
An editor who clicks the link in the manuscript's Code availability statement — or simply
asks for the DOI — finds nothing. That is a much worse position than arriving without a
DOI, because it reads as a false statement rather than an unfinished step.

**The parenthetical fallback does not rescue it**, and it introduces a second problem:
"The analysis code … *will be deposited in a public repository*" reads as "after
acceptance". Nature Portfolio's code policy requires code to be available **to editors and
reviewers**, not to readers after publication.

**Fix, in order of preference:**

**(a) Push the repository and publish the release before sending the letter.** Then the
present tense is true and the manuscript's Code availability statement resolves. This is a
single afternoon's work and is the recommended route:

```
GitHub: create a public repository, do NOT initialise it
git config user.name / user.email      # the current commit carries a placeholder author
git commit --amend --reset-author --no-edit
git remote add origin https://github.com/ACCOUNT/REPO.git
git push -u origin main
Zenodo: Settings > GitHub > toggle the repository ON
GitHub: Releases > Draft a new release > tag v1.0.0 > Publish
Zenodo: open the new record, copy the version DOI
```

**(b) If the letter must be sent first**, describe the state accurately:

> The analysis code is available to editors and reviewers at
> https://github.com/ACCOUNT/REPO and will be archived with a permanent DOI before
> publication; the Code availability statement gives the link.

Replace with the corrected paragraph below once the DOI exists.

**Corrected paragraph (use after (a)):**

> **Data and code availability.** All datasets are publicly available under the GEO
> accessions and the eQTL and GWAS resources listed in the Data availability statement.
> The analysis code has been deposited in a public repository and archived with a permanent
> DOI, both cited in the Code availability statement; the repository documents every
> analytical parameter, including the parameters of the robustness audit.

---

## ⚠️ 2. Two claims the repository contradicts

### 2a. "every analytical decision capable of biasing candidate selection having been fixed in advance"

**Current text (fourth paragraph):**

> Our reasons for selecting *Scientific Reports* are threefold. The analysis is technically
> sound and fully reproducible, **every analytical decision capable of biasing candidate
> selection having been fixed in advance.**

**The problem.** Two parameters were in fact adjusted during the analysis, and the
repository documents both:

- the robust rank aggregation score cutoff moved from **0.05 to 0.01**;
- the leave-one-dataset-out background of **20,000** was selected **by hand** after
  comparing three definitions (the six-cohort common-gene count of 11,982, per-cohort
  counts, and 20,000) — the original code comment reads `use_N20000 <- TRUE  # 校准后手工指定`.

Neither adjustment is misconduct, and both are recorded. But the claim as written invites
an editor to test it against the repository, and it does not survive that test. A
falsifiable overstatement about pre-specification is exactly the kind of thing that turns a
sound paper into a credibility question.

**Corrected sentence:**

> The analysis is technically sound and fully reproducible: every analytical parameter is
> declared in a single configuration script in the repository and tabulated in its
> documentation, and the robustness audit is deterministic, so its results regenerate
> exactly.

This is stronger, not weaker — it makes a claim the repository can substantiate.

### 2b. "including pre-specified parameters and random seeds"

**The problem.** Three separate inaccuracies in six words:

| Claim | Reality |
|---|---|
| "pre-specified parameters" | some were adjusted, see 2a |
| "random seeds" — plural | there is **one** seed in the project: `set.seed(20260914)`, for the GSE135251 AUC bootstrap |
| seeds at all, in this context | the thing the letter is implicitly pointing at — permutation testing — **does not exist** in this study |

The manuscript's Code availability statement is being edited to remove its permutation
clause for the same reason. If the letter keeps promising seeds, the letter and the
manuscript disagree, and a reviewer comparing them has to decide which one is wrong.

**Corrected clause:** see the replacement paragraph under item 1 — "the repository
documents every analytical parameter" makes the same point without the seed claim.

---

## ⚠️ 3. Two phrases in the summary that read as pre-registration

**Current text (second paragraph):**

> We integrated eight public bulk transcriptomic cohorts, six hepatic and two gastric, by
> robust rank aggregation, and subjected the shared candidates to a robustness audit
> **specified in advance**.

**Current text (third paragraph):**

> **First, all candidates were specified in advance** and were not re-selected after the
> discovery stage; validation was conducted in external cohorts that played no part in
> candidate selection.

**The problem.** "Specified in advance" and "specified in advance" in this position read as
pre-registration. There was none. The candidate set is an **output** of the discovery
screen: the robust rank aggregation with its score cutoff and direction-consistency filter,
intersected with the gastric differentially expressed genes. Nothing was fixed before the
data were seen.

The claim that **is** true, and that is the genuinely strong design feature, is *held-out
validation*: GSE135251 and GSE153224 played no part in candidate selection, and the
candidates were interrogated there once, in the pre-specified direction, without
re-screening. That is checkable from the code and is worth saying plainly.

**Corrected sentences:**

> We integrated eight public bulk transcriptomic cohorts, six hepatic and two gastric, by
> robust rank aggregation, and carried the shared candidates into a pre-specified
> robustness audit.

> **Two design features relevant to the interpretation of the work.** First, the candidate
> set was carried forward from the discovery stage and was not re-selected; validation was
> conducted in external cohorts that played no part in candidate selection.

"Pre-specified" survives in the second version only where it is accurate — the audit was
specified before it was run — and disappears where it was not.

---

## Everything else in the letter checks out

| Item | Status |
|---|---|
| Date — 29 September 2026 | consistent |
| Title — the 17-word version with "metabolic dysfunction-associated steatotic liver disease" | correct, matches the decided title |
| Article type | Article ✓ |
| Five retained genes listed correctly — CDHR2, ANXA4, IL32, CADM2, LGALS3 | ✓ |
| "eight public bulk transcriptomic cohorts, six hepatic and two gastric" | ✓ matches the discovery set |
| STROBE-MR checklist mentioned as accompanying | ✓ required by *Scientific Reports*; make sure it is actually uploaded as a separate file |
| Three suggested reviewers with emails | ✓ matches the pool; no co-authorship conflicts found |
| Excluded reviewers line | ✓ *"We have no request for excluded reviewers."* is the correct phrasing for the 3-name limit |
| Prior discussion with an Editorial Board Member | ✓ |
| No competing interests / ethics not required / funding grants | ✓ plausible for a re-analysis of public data; the ethics rationale in the Methods must match this wording |
| Corresponding author contact | ✓ email and telephone complete |

---

## Apply in this order

1. **Decide the deposition timeline** (item 1). Everything else is independent of it, but
   the letter's third paragraph depends on the answer.
2. Fix items 2 and 3 in the `.docx`.
3. Make the manuscript's three edits: HPA version and access date in Data availability, the
   GSE174478 accession, and the Methods/Table-S02/Code-availability permutation wording.
4. Push, release, take the DOI, fill both statements.
5. Send.

The wording in `cover_letter_SR_v2.md` is already corrected and can be copied from there.
