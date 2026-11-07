# Tm: two columns, and why two tools disagree

Which temperature to anneal to, which number matches which tool, and why a third
column was removed. [← back to the README](../README.md)

## Tm, and why two tools disagree

**Two** Tm are reported per primer, and each matches something exactly: `Tm PCR`, the
SantaLucia 1998 nearest-neighbour value the design was selected on, and `Tm UCSC`, the number
UCSC in-silico PCR prints for the same oligo. What matters first is the conditions, which are
stated rather than inherited from a library default:

| | default | flag |
|---|---|---|
| monovalent cations (Na⁺) | 50 mM | `--salt` |
| divalent cations (Mg²⁺) | 1.5 mM | `--salt` |
| total dNTP | 0.6 mM | `--salt` |
| annealing oligo | 50 nM | `--salt` |

These were previously left to `primer3-py`'s defaults, which is why the reported Tm is a
**PCR-condition** Tm. The Mg²⁺ term is worth about **6 °C** on a 20-mer: `PRKCG_E5F`
`ACTGTGGCTCCCTCCTCTAC` reads **60.3 °C** here and **54.4 °C** with Mg²⁺ and dNTP at zero, from
the very same tables. Most simple web calculators model Na⁺ alone, so they read lower — that is
a difference of convention, not a disagreement about the primer. Expect it when you paste a
primer into IDT or a similar tool, and do not anneal to what they tell you.

Matching the salt is not always enough, because the *tables* can differ too. UCSC in-silico PCR
prints a Tm from the older Breslauer 1986 parameters with the Schildkraut–Lifson salt term —
"the formula by Rychlik, Spencer and Rhoads 1990", i.e. Primer3's `PRIMER_TM_FORMULA=0` — while
Primer3 has defaulted to SantaLucia 1998 since 2.0. At identical salt the two disagree by
**5–8 °C on a 20-mer**, and the gap is sequence-dependent, so no constant offset converts one
into the other: for `POLG_E3` the forward primer reads 54.2 °C (SantaLucia) vs 59.3 °C (UCSC)
and the reverse 54.4 vs 62.0, which also inflates the apparent ΔTm of the pair from 0.0 to
2.7 °C. So `F Tm UCSC` / `R Tm UCSC` reproduces the UCSC number exactly, and cross-checking a
pair on hgPcr no longer looks like a bug. Data built before the column gains it from
`scripts/backfill_tm_ucsc.py`, with no redesign: the value depends only on the primer sequence.

A pair's **ΔTm** is stable across *salt* conventions but not across *tables*. `PRKCG_E5` is
0.2 °C apart under PCR conditions and the same 0.2 °C with the magnesium removed — salt moves
both primers together — but **5.6 °C** apart under UCSC's tables (58.3 vs 63.9), which is the
difference between a matched pair and one you would redesign. Simple GC rules fail in the other
direction: both primers of that pair are 60 % GC, so any GC-only formula calls them identical
whatever their sequences actually do. Matching is done on `Tm PCR`, and that is the column to
compare from pair to pair.

### Which temperature to use

**`Tm PCR`.** Always, for anything you do at the bench — it is `F Tm PCR (°C)` in the workbook
and the plain `Tm °C` column on the site. `Tm UCSC` is there for one purpose: so that a primer
pasted into hgPcr gives a number you recognise.

Each of the two matches something exactly, which is the rule that decided how many to report:

| column | tables | salt | matches, exactly |
|---|---|---|---|
| **`Tm PCR`** | SantaLucia 1998 | 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP | the value the design was selected on, under the conditions in the tube |
| `Tm UCSC` | Breslauer 1986 + Schildkraut | 50 mM Na⁺, no Mg²⁺ | what UCSC in-silico PCR prints |

`Tm PCR` wins on both axes at once: SantaLucia is the modern, better-validated set — which is
why Primer3 made it the default — and the PCR salt is what a standard reaction contains. It is
also, exactly, the number the design was chosen on: `PRIMER_OPT_TM`, the 57–63 °C window and
`PRIMER_PAIR_MAX_DIFF_TM` all act on it and on nothing else. So the usual rule of thumb — anneal
a few degrees below the lower primer's Tm, or touch down to it — applies to this column and to
no other number on the page.

`Tm UCSC` is a legacy display, and it is worth knowing it is *only* a display: isPcr matches by
sequence, not by Tm, and the command-line tool prints no Tm at all (we run it with `-out=bed`,
which has no such field). Use it to confirm you are looking at the same oligo UCSC is, never to
set a temperature.

**A third column used to sit between these two**, the same SantaLucia tables with the magnesium
removed (54.2 °C for the POLG forward primer). It was dropped, because it is the one number here
that matched nothing. It was not the design value — no primer was ever selected on it — and it
did not reconcile with the calculators it was meant for either: IDT defaults to 250 nM oligo,
not our 50 nM, and reports 56.5 °C for that same primer. So it agreed with neither us nor them,
while reading lowest of the three — which made it the number a hurried reader would anneal to,
5–6 °C too cold, with the non-specific products that follow. The magnesium gap it was explaining
is explained above in prose instead, where it cannot be mistaken for a target temperature.
`tm_na()` survives in the code, unpublished: designing under Na⁺-only conditions must land the
reported Tm exactly on it, which is the cleanest proof that the salt settings reach Primer3.

**Reproducing the UCSC column has a trap in it.** Switching `calc_tm` to the Breslauer tables
*without* also zeroing `dv_conc` and `dntp_conc` does not give UCSC's number: it gives UCSC's
number plus **8.56 °C**, on every primer. The offset is exact and
sequence-independent, because Primer3 converts divalents with `120·√(dv − dntp)` — no sequence
in it — and Schildkraut then adds a flat `16.6·log10([Na⁺])`:

```
16.6 · log10((50 + 120·√(1.5 − 0.6)) / 50) = 8.5566 °C
```

Measured over 22,654 distinct primers from the current build, min = max = +8.557. That is why
`tm_ucsc()` passes all four of `mv_conc`, `dv_conc`, `dntp_conc` and `dna_conc` explicitly
rather than relying on any default. The same arithmetic is why the design Tm cannot be converted
into UCSC's by adding a constant: SantaLucia's correction enters through ΔS scaled by
`0.368·(N−1)`, so it depends on length and sequence, and across those same 22,654 primers the
gap ranges from −7.8 to +5.5 °C.

The practical payoff of designing to one target is a single programme for the whole panel.
On the current build **97.9 %** of best pairs have both primers inside the 57–63 °C window
(median 59.8 °C) and 99.8 % are matched to within 2 °C, so one annealing temperature — 58–60 °C,
or a 65 → 58 °C touchdown — runs essentially every amplicon. The exceptions are the 5.3 % of
amplicons designed with relaxed settings, which are flagged **Check** precisely so they are not
loaded onto a plate on the assumption that they behave like the rest.

### Matching the two primers of a pair

`--max-tm-diff` (default **2 °C**) is the ceiling, and `--wt-tm-diff` (default **1**) is how hard
Primer3 prefers a matched pair. The weight matters more than the ceiling: **Primer3's own default
for it is 0**, so without it a matched pair is merely *allowed*, never *preferred* — which is why
10 % of the previous build's best pairs exceeded 2 °C.

Both levers change which pair wins, so they change primers. Measured on 179 amplicons across 10
genes, against the previous settings (ceiling 3, weight 0):

| settings | primers unchanged | amplicons with ΔTm > 2 °C | median ΔTm |
|---|---|---|---|
| ceiling 3, weight 0 (before) | — | 18 | 0.5 |
| ceiling 2, weight 0 | 91 % | 7 | 0.5 |
| ceiling 2, weight 1 (now) | 39 % | **0** | 0.2 |

The count of OK amplicons is the same in all three, so nothing is lost either way — it is a
choice between matched pairs and leaving a published panel alone. `--wt-tm-diff 0` keeps 91 % of
the primers and still removes most of the wide pairs, if that trade suits you better.

### On the site: alternatives re-ordered by ΔTm

Every amplicon view shows **ΔTm** for the pair, and the checkbox above the search box re-orders
the alternatives already published so that a pair matched to within 2 °C comes first. It changes
nothing in the data — Primer3's own rank 1 stays listed, and a pair is never promoted if it would
trade away genome-wide specificity. It matters much less than it used to: on the current build
only **187 of 86,437** best pairs (0.2 %) exceed 2 °C, against 10 % before `--wt-tm-diff`, and
141 of those 187 have a matched alternative behind the checkbox.

