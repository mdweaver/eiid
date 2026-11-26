# Proof review: PHASE6a (M-COUNT, build step 1; register rows R19, R27)

*Independent proof review of `docs/extension/PHASE6a_count_layer.md` (cited as 6a) at branch `claude/mcount-theory`, head `d219d436`. Proofs only; the R code is reviewed separately. The theory documents are not edited here. Checked against `notation.tex`, `PHASE2c_monotone.md` (Def. 1, 2; M1–M3; Remark M2.1), `PHASE3_inference.md` (§0.1–0.2, Thm P3.5(a), (b2), Prop. P3.5(c1), Thm P3.8), `PACKAGE_DESIGN.md` §8 and `next_steps/NEXT_STEPS_RECOMMENDATION.md` (NS) step 1 and Part 3A.*

*Reviewer's checks: `phase6a_review_checks.py` (output `phase6a_review_checks_output.txt`, 13/13 PASS, ~25 s). Tags there: **[CE]** a counterexample to a statement as written in 6a; **[V]** an independent verification (different generator from the author's); **[B]** a boundary case showing that a hypothesis of 6a is needed (no error in 6a).*

*Reproduction of the author's script: `verify_phase6a.py` rerun with Python 3.13, scipy 1.18.1, numpy 2.5.3: output **identical** to `verify_phase6a_output.txt` (32/32 PASS, ~45 s).*

---

## 1. Verdict in brief

The core of 6a is sound. The count-layer LP (R19(a)), Duncan–Davis as its margins-only case (Lemma 2.1), the three TU classes (R19(b), conclusion), the MILP fallback (R19(c)), the witnesses (R19(d)), non-integer margins (R19(e)), the count form of the realised associations (Lemma C1) and the cross-layer coverage and test (R27(i)–(ii)) are all correct in substance. I checked every Ghouila-Houri partition line by line, and the integrality conclusions hold in 1,800 further random designs with deeper laminar families and infeasible right-hand sides (RC2).

Six things must change before the register rows can lose "proof review pending":
* one false sentence in the statement of R19(b) (issue 1);
* one false sentence and one missing scope hypothesis in Prop. C2 (issues 5, 6);
* a missing probabilistic formulation of $H_C$ in R27 (issue 2, the most substantive);
* a missing measurability hypothesis in Lemma 4.1 (issue 3);
* an unproved row of the §4.3 status table (issue 4).

Each has a short fix that leaves the results' use in the package unchanged, except issue 2, which changes what `print()` may say when a belief row about the realised table is in the count layer.

| Claim | Verdict | Issues |
|---|---|---|
| Lemma 1.1 (rounded system) | accept | 12 (contract) |
| R19(a) (LP, valid, sharp for continuous tables) | accept | 14 |
| Lemma 2.1 (DD = margins-only LP for single cells, one-column contrasts; cellwise outer otherwise) | accept | 9 |
| R19(b) (TU classes (i)–(iii); integer endpoints = rounded LP; intermediate integers) | **accept with fixes** | **1**, 8 |
| Corollary 2.2 | accept with fixes | 7, 8 |
| R19(c) (outside the classes: LP outer, MILP exact; strict gaps occur) | accept with fixes | 7, 15 |
| R19(d) (witnesses) | accept | — |
| R19(e) (non-integer margins) | accept | 12 (contract) |
| Lemma C1 (associations in counts) | accept | 14 |
| Prop. C2 (2c M2/M3 in counts; integer lattice; integer refutation) | **accept with fixes** | **5**, **6**, 13 |
| Lemma 4.1 (residual allowance $S_\varphi$) | **accept with fixes** | **3** |
| R27(i), (ii) (coverage; disjointness test) | **accept with fixes** | **2**, **4** |
| R27(iii-a) (not an identified set) | accept with fixes | 10 |
| R27(iii-b) (not minimal; computed example) | accept as a computed example | 11 |
| §9 claimed errors in other documents | all real or accurately described | 9 (extends item 3) |

Severity: **error** = a statement that is false as written; **gap** = a missing hypothesis or a missing step, where the claim is true once it is added; **wording** = imprecise or incomplete, with no false content.

---

## 2. Numbered issues

### Issue 1 (error, statement of R19(b); the conclusion stands). "The constraint matrix of the rounded system is TU" is false in class (iii).

The theorem states, for all three classes, "the constraint matrix of the rounded system (margins and totals; bounds are not rows) is TU". The proof of class (iii) (step 5) does **not** show this. It eliminates $n_{i2c}=b_{ic}-y_{ic}$ and shows that the **reduced** system in $y$ is TU (unit rows form one laminar family, the totals of all columns another). That proves that $F_{\rm rd}$ is an integral polytope, which is all that items 1–3 need. But the original margins-plus-totals matrix need not be TU.

*Counterexample (RC1).* Three 2×2 units. In column 1, totals on cell (1,1) over $\{0,1,2\}$ and on cell (2,1) over $\{0,1\}$; in column 2, a total on (1,2) over $\{1,2\}$. Each column's sets form one laminar family, so this is class (iii). The rows $T_{\{012\},(1,1)}$, $E^{\rm row}_{2,1}$, $T_{\{12\},(1,2)}$, $E^{\rm col}_{1,2}$, $E^{\rm row}_{1,2}$, $T_{\{01\},(2,1)}$, $E^{\rm col}_{0,1}$ and the variables $n_{0,11},n_{2,11},n_{2,12},n_{1,12},n_{1,22},n_{1,21},n_{0,21}$ form an odd cycle: the 7×7 minor is 2. On 300 random right-hand sides for this design the rounded LP nevertheless equals the MILP and every vertex is integral, as the substitution argument says.

*Why the GH signing of steps 1–4 cannot be reused.* With totals in several columns, the per-unit orientation $\varepsilon_i$ must serve every column at once. When $E^{\rm col}_{ic}$ is chosen and $E^{\rm col}_{ic'}$ is not, no single $\varepsilon_i$ keeps both columns' variables in $\{-1,0,1\}$ once the totals add $\tau\in\{0,1\}$.

*Fix.* Replace "Then the constraint matrix of the rounded system … is TU, and $F_{\rm rd}$ is an integral polytope" by: "Then $F_{\rm rd}$ is an integral polytope. In classes (i) and (ii) the margins-plus-totals matrix is TU. In class (iii), after the integral affine substitution $n_{i2c}=b_{ic}-y_{ic}$, the matrix of the reduced system is TU (the original matrix need not be; RC1)." In the proof's first paragraph, "it suffices to show that the matrix of $E$ and $T$ rows is TU" should be restricted to (i)–(ii), with step 5 named as the argument for (iii). Item 3's step 6 ("adding the row keeps the matrix in the class") is fine for (iii) when read in the reduced system. The register row R19(b) says "TU for …"; it should say "integral (TU, or TU after eliminating one row when $R=2$)".

### Issue 2 (gap, R27(i)–(ii)). What "under $H_C\wedge H_M$" means when a count row is a belief about the realised table.

Under $H_M$ (P3-S1a) the realised table $n^\circ$ is **random** given the design; 6a's own (iii-a) says so. A count row that is a **data identity** (the margins, a published total, an audited known cell: rows whose right-hand side is the realised value itself) holds in every replication. A count row that is a **belief** (a realised ordering, a realised sign, a box, a stated floor on a count) is an event $\{n^\circ\in F\}$ whose probability under $H_M$ is generally strictly between 0 and 1. Step 1 of the proof ("under $H_C$, $B_\varphi\in\mathcal C_\varphi$ surely, so the events coincide") uses $H_C$ as a sure statement. Two readings are then possible, and neither gives the theorem as written:
* *Conditional on $H_C$:* $P(B_\varphi\in\mathcal C_\varphi\cap\mathrm{CI}_\varphi\mid n^\circ\in F)=P(B_\varphi\in\mathrm{CI}_\varphi\mid n^\circ\in F)$, which can be far below $1-\alpha$. Likewise the conditional rejection rate can be far above $\alpha$. RC8: pooled 2×2 LCE true, an exact oracle 95% interval for $B$, and a belief row "$B\ge t$". With $t=\tilde B+1.0\,\mathrm{sd}(E)$, $P(\text{row})=0.161$ and the conditional coverage of $\mathcal C\cap\mathrm{CI}$ is 0.838. With $t=\tilde B+2.4\,\mathrm{sd}(E)$, $P(\text{row})=0.008$ and the conditional rejection rate is 1.000.
* *Unconditional:* then $H_C$ is not a hypothesis that "holds"; it fails with probability $P_{H_M}(n^\circ\notin F)$.

*What is true (same one-line proof).* Under $H_M$ alone:
$$P(B_\varphi\notin\mathcal C_\varphi\cap\mathrm{CI}_\varphi)\le P(B_\varphi\notin\mathrm{CI}_\varphi)+P(n^\circ\notin F),\qquad P(\text{reject}\cap\{n^\circ\in F\})\le P(B_\varphi\notin\mathrm{CI}_\varphi).$$
* When every count row is a data identity, or more generally $P_{H_M}(n^\circ\in F)=1$, the second term vanishes and (i)–(ii) hold exactly as stated. This covers the margins-only case, published totals and audited known cells, which are the cases X2 simulates.
* With belief rows, (i) holds with level $1-\alpha-P(n^\circ\notin F)$, which is unknown. (ii) bounds the probability of "rows true and rejected" by α; it is not a test with conditional level α given the rows. This parallels Remark 4.2's union bound for estimated margins, with $P(n^\circ\notin F)$ in the role of $\alpha_C$.

*Fix.*
1. Define $H_C$ in §4.1 as "$P_{H_M}(n^\circ\in F)=1$", or split the count rows into identity rows and belief rows.
2. Restate (i)–(ii) with the two displayed inequalities, noting that they reduce to the present statements for identity rows.
3. Change the second bullet of (ii) to: "with belief rows, it bounds the probability of rejecting while the rows are true; a rejection is evidence against $H_M$ only with the rows maintained, and its level is α only if the rows are taken as certain".
4. Contract C8 must print that qualification when a belief row is present. The current "rejects (model + its rows) at level α given the count rows" overstates it, because "given" reads as conditioning.
5. §9 item 2 (NS Part 3A) remains correct as a criticism of NS; the suggested replacement wording "given the count rows" should be "with the count rows maintained".

The register row R27 must carry the same restriction.

### Issue 3 (gap, Lemma 4.1). The coefficients must be design-measurable, and $S_\varphi$, $E_\varphi$, $\tilde B_\varphi$ depend on how the target is written.

Step 1 writes $E_\varphi=\sum_i(Z_i-E[Z_i\mid V_i])$ with independent $Z_i$. This needs $v_{irc}$ to be fixed given the design $\{(N_i,X_i,D_i)\}$ (P3-S2). The outcome margins $b_{ic}=N_iT_{ic}$ are **random** under P3. Two consequences:
* *(a) Targets whose coefficients involve outcome margins are not covered.* An example is a share among the $c$-voters, $Y_{rc}/Y_c$, with $v=1/Y_c$. Then the $Z_i$ are not independent sums and $E[E_\varphi\mid\text{design}]\neq0$ in general. §1.2 allows "any linear functional with known coefficients" and C2 allows "any linear count functional over a known denominator"; the lemma needs "known given the design".
* *(b) Representation dependence (RC7).* An intercept $\varphi_0$ may use observed outcome margins: it is observed, and the P3.5 interval for $\sum v a m$ simply shifts by it. But then the split $B_\varphi=\tilde B_\varphi+E_\varphi$ depends on the representation. 6a's own form $D^r_c=\kappa_rY_{rc}-Y_c/(N-A_r)$ and the cell form $Y_{rc}/A_r-\sum_{s\ne r}Y_{sc}/(N-A_r)$ give different $E_\varphi$: they differ by $(Y_c-E[Y_c\mid\text{design}])/(N-A_r)$. They also give different $S_\varphi$: 0.0465 vs 0.0555 in RC7's design. Each is valid with its own allowance. The variance bound itself is fine: RC7 checks the exact variances.

*Fix.* Add to Lemma 4.1: "$v_{irc}$ are functions of the design only; $\varphi_0$ may depend on observed data." Then add a remark: "$\tilde B_\varphi$ and $S_\varphi$ refer to the stated representation. Adding $\lambda_{ic}(\sum_rn_{irc}-b_{ic})$ with design-measurable $\lambda$ gives another valid representation, so the package may use any fixed one, e.g. the one minimising $S_\varphi$, chosen from the design alone. Targets with outcome-dependent denominators are not covered." §5 item 1 ("$\tilde B_\varphi$ is not a function of the table") should say "for a fixed representation".

### Issue 4 (gap, §4.3 table and register row R27). The P3.8 row is proved only for single cells.

P3.8's finite-sample residual allowance is $t_E=(\tfrac12\log(2/\alpha_E)\sum_iw_{ir}^2)^{1/2}$, by Hoeffding for $B_{rc}$. Lemma 4.1 extends only the CLT/Popoviciu allowance. For multi-cell targets the P3.8 row therefore has no proof in 6a. The extension is one line: $Z_i$ lies in an interval of length $\rho_i$ (Lemma 4.1 step 2) and the $Z_i$ are independent, so Hoeffding gives $t_E(\varphi)=(\tfrac12\log(2/\alpha_E)\sum_i\rho_i^2)^{1/2}$. The estimation part of P3.8 already holds for any linear functional.

*Fix.* Add this as Lemma 4.1(b) (with issue 3's measurability hypothesis), or restrict the P3.8 row to single cells.

### Issue 5 (error, Prop. C2(b), one sentence). The integer gap set is not "the lattice inside the continuous interval".

C2(b) says that the integer values of $Y_{rc}$ are all integers in $[\max(Y^-_{\rm rd},\lceil\cdot\rceil),\min(Y^+_{\rm rd},\lfloor\cdot\rfloor)]$. That is correct. It then says: "So the integer gap set is the lattice $\{D^r_c(k)\}$ … inside the continuous interval." This is false when fractional boxes sit on several units.

*Counterexample (RC5), the T5 design.* Two units ($N=5$, $V=1$, yes $=1$) with the realised ordering boxes $n_i\ge0.2$ and a between sign $Y\ge Y^{NM}=0.4$. The continuous range of $Y$ is $[0.4,2]$, which contains the integers $\{1,2\}$. The integer values are $\{2\}$. The formula with $Y^-_{\rm rd}=2$ is right; the sentence is not. The sentence is true only without unit rows, or with integer boxes, which is 2c Remark M2.1's case.

*Fix.* "So the integer gap set is $\{D^r_c(k)\}$ for the integers $k$ in that interval. Without unit rows, or with integer boxes, this is the lattice inside the continuous interval; with fractional boxes on several units it can be strictly smaller (T5)."

### Issue 6 (gap, Prop. C2(a)–(c) scope). One sign set; the target is the aggregate over that set; CR without slack.

* C2's thresholds $Y^{NM}_{rc}$ and $Y^{ER}_{rc}$ are computed over **one** unit set $\mathcal S$ and act on $Y_{rc}(\mathcal S)$. C2(a) ("any system of one-vs-rest sign rows … cuts $Y_{rc}$ to an interval $J$") and C2(b) (the formula for "the integer values of $Y_{rc}$") are correct only when every sign row is over the same $\mathcal S$ and the target is $Y_{rc}(\mathcal S)$. C2(b) allows sign sets forming a laminar family, and its proof uses "$Y_{rc}$ over all units". With sign rows over several sets, or a target over a different set, the formula is wrong. RC5b: a between sign over $\mathcal S=\{0\}$ gives the lower end 1 for $Y_{11}$(all units), while $Y^{NM}$ over all units is 2. The integer statements of R19(b) items 1–2 (rounded LP = integer endpoints) remain true for such systems; only the closed form and the lattice description need the single set.
* C2(a) says "CR is $J_+\cup J_-$, one of which is a point or empty". This holds without slack (2c Def. 2's CR has none). With slack dials on the CR rows, both branches can be nondegenerate intervals (for $Y^{NM}=Y^{ER}$ and $\varepsilon>0$).

*Fix.* Add to C2's preamble: "all sign rows are over one unit set $\mathcal S$; $Y_{rc}$ means $Y_{rc}(\mathcal S)$". For several sign sets, refer to R19(b) items 1–2 (integer endpoints by the rounded LP) without the closed form. Add "(without slack)" to the CR clause.

### Issue 7 (wording, R19(c), §3.4, register row R19). "Totals on cells in different rows and columns ($R\ge3$)" needs $R\ge3$ **and** $C\ge3$.

With $C=2$ the transposed class (iii) applies: for each row, the sets on that row's two cells form one laminar family. RC10: 1,500 designs of 6a's own C1 type with $C=2$ gave 0 fractional or LP≠MILP endpoints. 6a's C1 example is 3×3, so nothing computed is wrong; only the qualifier is. The same applies to "one-vs-rest signs or totals on cells in different rows and different columns over overlapping units ($R\ge3$)" in §3.4.

*Fix.* "$R\ge3$ and $C\ge3$" in all three places.

### Issue 8 (wording, Corollary 2.2, §3.4, contract C4). Class (ii) for signs, and the classifier must test one class for the whole row set.

* Corollary 2.2 lists one-vs-rest signs on one cell (class (i)) and, for $R=2$, several columns (class (iii)). It omits what class (ii) gives directly: one-vs-rest signs and totals on **several rows of one column**, or several columns of one row, with one laminar family per cell. §3.4 implies this only by contrast.
* R19(b) requires every non-union row to fall in **one** class. Mixtures are not covered and can fail. RC4: two families on one cell plus a total on another cell of the same column gave LP −7 vs integer −3; $R=2$ with crossing sets on one column was LP-feasible and integer-infeasible. 6a states the classes correctly, but contract C4 ("if every row is in the classes of R19(b)") reads as if each row could be checked separately.

*Fix.* Add the class (ii) bullet to Corollary 2.2. In C4, write "if the whole set of non-union rows lies in one class (per CR branch)". Note for the code step that "two laminar families" is decided by 2-colouring the crossing graph of the sets (bipartiteness).

### Issue 9 (wording, extends §9 item 3 and contract C2). Name the one-vs-rest gap $D^r_c$ explicitly.

For $R\ge3$, $D^r_c=B_{rc}-\sum_{s\ne r}\tfrac{A_s}{N-A_r}B_{sc}$, written over cells, is a "several same-sign cells of one column" target, so the package's cellwise `duncan_davis` is outer for it. Its exact margins-only range is the single-cell DD range of $Y_{rc}$ mapped by $D=\kappa_rY_{rc}-Y_c/(N-A_r)$ (Lemma C1(a) with Lemma 2.1(a)). RC6: strictly wider in all 200 random designs, by up to 0.87. Since $D^r_c$ is the M-SIGN default contrast, this is the case users will hit first.

*Fix.* Add it to §9 item 3 and to C2 ("must equal the current `dd_lower`/`dd_upper` for single cells and two-row one-column contrasts; for $D^r_c$ with $R\ge3$ it is narrower").

### Issue 10 (wording, R27(iii-a)). "Identified set".

The argument shows that no table of $F$ is excluded by $H_C\wedge H_M$: the two-point law has mean $m(V_i)$ and puts mass on the table's rates. It needs $m(V_i)$ in the relative interior for **every** unit (or the stated structural-zero variant), with $m$ evaluated at some $k$ in the model set. That is a statement about the values of the realised $B_\varphi$ compatible with the observables and the hypotheses, not an identified set in the usual (population-law) sense. With a single realisation the latter is not well defined here.

*Fix.* "The set of values of $B_\varphi$ compatible with the observed margins and $H_C\wedge H_M$ is $\mathcal C_\varphi$ (when every $m(V_i)$ is interior): $H_M$ excludes no table of $F$." Keep "(iii-a) P" with that wording.

### Issue 11 (gap, minor, Remark 4.4 / (iii-b)). The set $J$ is valid only if the known-cell units are selected by design.

$\mathrm{CI}_{\mathcal U}$ is a model interval over the units without known cells. Its coverage needs $\mathcal K$ to be chosen from the design only, as in P3's design-based selection result (Prop. P3.7(b)), plus the Lindeberg condition over $\mathcal U$. Audits targeted on outcomes (close or anomalous precincts) violate this. X5 draws $\mathcal K$ at random, so the example stands.

*Fix.* State the hypothesis in Remark 4.4. (iii-b) remains a computed example (see §4 below).

### Issue 12 (wording, contract C4 and C7; Lemma 1.1, R19(e)). Exact arithmetic for rounding and for "integer margins".

Lemma 1.1 is exact, but its use needs exact $\lceil\cdot\rceil$/$\lfloor\cdot\rfloor$. RC9: in 4,008 random small populations where $Y^{ER}$ is an integer in rational arithmetic, the double-precision floor was wrong in 1,281. A wrong floor on an upper threshold removes an attainable integer, so "integer-exact" would be false. Similarly, "auto" must decide whether margins are integers, and the package's reference data use $N\cdot X$, which can be near-integers in floating point.

*Fix.* C4: "thresholds are computed in exact rational arithmetic (or rounded with a certified tolerance, with the tolerance printed); margins count as integers only if they are exactly integer in the input".

### Issue 13 (wording, Prop. C2(c)).

Without unit rows, $Y^\pm$ are integers, so "$\lfloor Y^{ER}\rfloor<Y^-$" is equivalent to "$Y^{ER}<Y^-$". A single sign is integer-refuted iff it is continuously refuted, and only CR (or several signs together) can be integer-infeasible while continuously feasible. Saying so would make the content clearer. No error.

### Issue 14 (wording, §1.1, §1.3, Lemma C1). Degenerate units.

Units with $N_i=0$ should be dropped explicitly: Lemma C1 divides by $N_i$. The odds-ratio dial needs a convention when a margin is 0 or the ratio is $0/0$ or $\infty$. The ordering box needs $0<a_{ir}<N_i$ (stated). $\delta^r_W$ needs $\operatorname{Var}_p(X_r)>0$ over the **sign set** $\mathcal S$ (stated for C2, implicit in Corollary 2.2).

### Issue 15 (wording, R19(c) examples).

The pairwise $\delta_B(r,s)$ example (LP 45.098, integer 45) shows a fractional vertex of $F_{\rm rd}$, but its integer endpoint equals the floor of the LP endpoint. R19(c)(3)'s "integer endpoints beyond the rounded LP endpoints" rests on the pairwise-ordering, two-cell and crossing-totals examples and on N10c. Say so, or replace the $\delta_B(r,s)$ instance with one whose integer endpoint differs from the rounded one (C1's search can be pointed at that).

---

## 3. Line-by-line notes on the main proofs (no issue unless numbered above)

* **Lemma 1.1.** Chvátal rounding of one row preserves the integer points; the inward rounding implies $F_{\rm rd}\subseteq F$. Correct.
* **R19(a).** Steps 1–5 correct. Item 4's slack statement is 2c M2(d) verbatim.
* **Lemma 2.1.**
  * Step 2's projection argument ([NW] on the residual $R\times(C-1)$ problem) is correct.
  * Step 3: both remainders are $\le N_i-a_{ir}$. Checked.
  * Step 4: I rechecked both cases. If $a_{ir}\le b_{ic}$, the conditional minimum is $b_{ic}-N_i+a_{is}$. If $a_{ir}>b_{ic}$, it is 0, and the Fréchet minimum is also 0 because $b_{ic}+a_{is}<a_{ir}+a_{is}\le N_i$. Correct.
  * Step 5's example: reproduced by hand ($[0,0.5]$ vs $[0,1]$; $[0,1]$ vs $[0,2]$).
* **R19(b), GH partitions.**
  * Step 1: margins with a per-unit orientation give $\varepsilon_i(e^{\rm row}-e^{\rm col})\in\{-1,0,1\}$.
  * Step 2: depth signing of a laminar family gives an alternating sum along the chain of sets containing $i$, so $\{0,1\}$. Correct; any strictly containing chosen set lies in the chain.
  * Step 3, class (i): with $\tau_i\ne0\ne\delta_i$, $\varepsilon_i=-\tau_i\delta_i$ gives $\tau_i-\tau_i\delta_i^2=0$. The other variables of unit $i$ get margin terms only, and the choices are per unit. Correct.
  * Step 4, class (ii): both cases give $\{-1,0\}+\{0,1\}$ on column $c_0$. The other columns get margin terms in $\{-1,0,1\}$. Correct, and it needs exactly one family per cell (with two families $\tau\in\{-1,0,1\}$ and case 1 can reach −2; consistent with RC4).
  * Step 5, class (iii): the reduced system is correct. A total on $(2,c)$ maps to one on $y_{\cdot c}$ with integer right-hand side after rounding. The unit rows and the totals are two laminar families on disjoint-by-column supports. Correct for the reduced system (issue 1 concerns only the statement).
  * Boxes and known cells enter as bounds, rounded to integers, so [HK] applies.
  * Item 1 needs boundedness for "nonempty ⇒ integer vertex", which holds ($0\le n\le\min(a,b)$).
  * Item 3: a single unit and the set of all units are laminar with any family, so the examples of "addable" are right. RC3 shows that "addable" cannot be dropped: a crossing target set has integer values $\{0,2\}$.
* **R19(c)–(e).** Correct. RC4 and the author's C1/N10c confirm that strictness occurs.
* **Lemma C1.**
  * (b): matches 2c Def. 1.
  * (c): matches 2c M1(b). I checked $NM_r=Y^{NM}_{rc}/A_r$ and $\delta^r_B=\mu_r(B_r-NM_r)$ by hand.
  * (d): matches 2c M1(c). The threshold algebra $D\le ER+\varepsilon/\operatorname{Var}$ ⇔ $Y\le Y^{ER}(\varepsilon)$ is correct.
  * (e): the unit coefficient of $\delta^r_B$ in counts is $p_iX_{ir}(1-X_{ir})(1/a_{ir}+1/(N_i-a_{ir}))=1/N$, constant across units. That is the precise sense in which NS's "weight units differently" (N10c comment) is wrong.
* **Prop. C2.**
  * (a): the projection of the R×C problem onto $Y_{rc}$ equals that of the collapsed problem (both are summed Fréchet bounds).
  * The equal-rates table $n_{isc}=a_{is}b_{ic}/N_i$ lies in $\mathcal P$ and hits $Y^{NM}$. Correct (apart from issue 6).
  * (b)–(c): see issues 5, 6, 13.
* **Lemma 4.1.**
  * Steps 2–3 (simplex range, Popoviciu) are correct, including dependence within a unit.
  * Step 4 transfers P3.5(a) step 5's subsequence argument verbatim, using $\lvert Z_i-EZ_i\rvert\le\rho_i$.
  * $z_{\alpha_E}$ is two-sided as in P3.
  * See issues 3 and 4.
* **R27.** (i)–(ii) are correct for identity rows (issue 2). Remarks 4.1–4.3 are correct. Remark 4.3's "$\mathcal C_\varphi\subseteq$ DD" holds because the margins are rows of $F$.

---

## 4. The author's flagged questions

1. **General dominance for R27(iii-b).** It is not needed and should not be claimed. The statement as worded ("not in general the smallest valid set") is an existence claim. It is established by one design in which a valid set ($J\cap\mathcal C$; valid by P3.5 on $\mathcal U$, given issue 11's design-based selection) is narrower in the computed realisations. Keep the status "computed example", not P. Note that "narrower" is in width: $J\cap\mathcal C\not\subseteq\mathcal C\cap\mathrm{CI}$ in general.
2. **Integer theory outside the TU classes.** R19(c) claims only "LP outer, MILP exact", and that is proved; nothing more is needed for the contract. The class boundaries are genuine (RC4), and the corrected qualifier is issue 7.
3. **Bhatia–Davis for multi-cell targets.** It is correctly not claimed. Contract C8 must forbid it for multi-cell targets, not only prefer $S_\varphi$. §9 item 7's deferral of the package's current multi-cell allowance to the code review is right.
4. **X2's size 0.000.** It is uninformative, as the author says, because the interval is conservative. The proof carries (ii); the simulation adds nothing on tightness. If an illustration that the bound can be approached is wanted, use an exact oracle interval with identity rows only (RC8's set-up without the belief row).
5. **May "P" stand before review?** Under PD §8's standard, "**P**; proof review pending" is the right label for the period of review, and the register rows say that. After this review:
   * R19 may drop "pending" once issues 1, 5, 6, 7 and 8 are fixed in 6a (5–6 are in Prop. C2, which the R19 row cites).
   * R27 may drop "pending" once issues 2, 3, 4, 10 and 11 are fixed, with the register row restricted per issue 2.

---

## 5. §9 (errors claimed in other documents): verification

1. **NS §3.1 / Part 3A / N10c comment: real.**
   * NS line 202 ("Pooled realised signs and CR couple cells across units, so total unimodularity fails: realised CR is two MILPs") is false for the W1 one-vs-rest signs on one cell, by Lemma C1(e) and class (i). The author's S2 and my RC2 (class (i)) confirm it.
   * The N10c comment in `verify_next_steps.py` (rows "as the realised δ_B and δ_W rows of 2c M1 weight units differently") is wrong in count coordinates: see the coefficient computation in §3 above.
   * Caveat to add: NS's statement *is* true for pairwise/W2 signs and for signs on cells in different rows and columns with $R,C\ge3$ (issue 7).
2. **NS Part 3A first bullet ("does not reject the ordering"): real.** With the ordering in the count layer, disjointness involves it. See issue 2 for the right replacement wording.
3. **Package `duncan_davis` for combination targets: real.** `.rxc_interval_lin` (`rxc_infer.R` l. 141) is cellwise interval arithmetic. Extend with issue 9 ($D^r_c$).
4. **2c Remark M2.1 "Open (R15)": an accurate scope update**, not an error. The phrase "for $R=2$ with joint columns, the gap is exactly the rounding of the rows" is right in the sense of R19(b) item 2 (endpoints), not as a lattice description (issue 5).
5. **NS R19(b) wording: real.** TU is a property of the matrix, and integer endpoints need integer right-hand sides. The author's T5 confirms it.
6. **PD §4.2 footnote f:** no error, agreed.
7. **P3.5(a):** no error, agreed. Issue 4 adds the P3.8 analogue.

---

## 6. Changes required before the register rows move from "proof review pending"

**R19** (6a §§1–3, PD §8 row R19):
* [ ] Issue 1: restate R19(b)'s conclusion as integrality; TU of the original matrix only for (i)–(ii); (iii) via the substitution. Update the register wording.
* [ ] Issue 5: correct the lattice sentence in Prop. C2(b).
* [ ] Issue 6: one sign set $\mathcal S$ and target $Y_{rc}(\mathcal S)$ in Prop. C2; CR "without slack".
* [ ] Issue 7: "$R\ge3$ and $C\ge3$" in R19(c), §3.4 and the register row.
* [ ] Issue 8: class (ii) bullet in Corollary 2.2; "one class for the whole row set" in C4.
* Recommended: issues 9, 12, 13, 14, 15.

**R27** (6a §4, PD §8 row R27):
* [ ] Issue 2: define $H_C$; restate (i)–(ii) with $+P(n^\circ\notin F)$; distinguish identity rows from belief rows in the statement, in C8's printed text and in the register row.
* [ ] Issue 3: design-measurable $v$ in Lemma 4.1; the representation remark.
* [ ] Issue 4: Hoeffding allowance $t_E(\varphi)$ for P3.8, or restrict the §4.3 P3.8 row (and the register's "P3.8 P") to single cells.
* [ ] Issue 10: "compatible values" wording in (iii-a).
* [ ] Issue 11: design-based selection hypothesis in Remark 4.4.

None of these changes the LP the package solves, the TU routing (apart from issue 8's whole-set test), or the margins-only cross-layer test. Issue 2 changes what may be printed when belief rows enter the count layer. Issue 3 constrains which targets get an actual-rate model interval.

---

## 7. Round 2 confirmation (6a at `claude/mcount-theory` ed709f11)

*Confirmation pass on the author's fixes (6a §10). I merged ed709f11 into this branch as a fast-forward, read the full diff of 6a and PD §8, and reran both scripts.*

*Reruns:*
* `verify_phase6a.py`: **36/36 PASS**, output identical to the committed `verify_phase6a_output.txt`. The new checks T6, S5, X6 and G2 reproduce RC1, RC5 and RC8 and test Lemma 4.1(b).
* `phase6a_review_checks.py`: **14/14 PASS**. That is the 13 checks of round 1, unchanged and all still passing, plus a new check RC11.

### 7.1 Verdict per issue

| Issue | Verdict | Notes |
|---|---|---|
| 1 class (iii) matrix not TU | **resolved** | The statement now concludes integrality. TU of the original matrix is claimed for (i)–(ii) only, and (iii) is proved via the reduced system; the proof's first paragraph, step 5 and step 6 all say so. The register row is restated to match. T6 reproduces the minor 2. |
| 2 belief rows are events | **resolved** | See §7.2. |
| 3 design-measurable $v$; representation | **resolved** | The hypothesis is now in Lemma 4.1 and §1.2, and outcome-margin denominators are excluded. Remark 4.1a covers the representation question, requiring one fixed from the design alone. C2 and C8 are updated. |
| 4 P3.8 for multi-cell targets | **resolved** | Lemma 4.1(b) adds the Hoeffding allowance and the §4.3 row is restated; see §7.2. |
| 5 C2(b) lattice sentence | **resolved, with a residual one-phrase error (N1 below)** | The false sentence is replaced and S5 reproduces RC5. However, the replacement condition "without unit rows, or with integer boxes" is still too weak. |
| 6 one sign set; CR slack | **resolved** | The preamble fixes $\mathcal S$ and the target $Y_{rc}(\mathcal S)$. With several sign sets the note now falls back to R19(b) items 1–2. The CR clause says "without slack", and (b) requires $\mathcal S$ to be addable. |
| 7 $R\ge3$ and $C\ge3$ | **resolved** | Fixed in the R19(c) example, §3.4, §9 item 1 and the register row. |
| 8 class (ii) bullet; one class for the whole row set | **resolved** | Fixed in the R19(b) statement, Corollary 2.2 (bullet and proviso), C4 (the crossing-graph bipartiteness test is stated correctly) and C10. |
| 9 $D^r_c$ | **resolved** | Fixed in Lemma 2.1(c), C2 and §9 item 3. |
| 10 "identified set" | **resolved** | (iii-a) is restated as the set of compatible realised values, with interiority required for every unit at some $k$ in the model set. |
| 11 design-based $\mathcal K$ | **resolved** | Remark 4.4 now carries the hypothesis, citing P3.7(b) and Lindeberg over $\mathcal U$. |
| 12 exact arithmetic | **resolved** | C4 now requires exact rationals or a printed tolerance, and counts margins as integers only if they are exactly integer in the input. C10 tests it. |
| 13 C2(c) | **resolved** | |
| 14 degenerate units | **resolved** | $N_i>0$ is stated. The odds-ratio convention is correct: with all margins positive, the ratio is 0 at the lower end of the cell's range and ∞ at the upper end, never $0/0$. The variance condition is required over each sign set. |
| 15 $\delta_B(r,s)$ example | **resolved** | The example is restated as a fractional vertex only, and R19(c)(3) is reworded. |

### 7.2 What the fixes introduced

**Lemma 4.1(b) (Hoeffding): correct.**
* Given the design, the $Z_i$ are independent and each lies in an interval of length $\rho_i$; step 2 shows this holds whatever the dependence inside a unit, provided $v$ depends on the design only.
* Hoeffding for $\sum_i(Z_i-EZ_i)$ then gives $2\exp(-2t^2/\sum\rho_i^2)$, which equals $\alpha_E$ at $t_E(\varphi)$. For a single cell, $\rho_i=w_{ir}$, which is P3.8 step 6 exactly.
* $\varphi_0$ does not enter $E_\varphi$, so letting it depend on the data is harmless. The interval for $\tilde B_\varphi$ is $\varphi_0$ plus P3.8's interval for the design-coefficient functional $\sum v a m$.
* G2 checks the extremal law.

**R27 union bounds $(*)$: correct.**
* Step 1 holds pointwise on $\{n^\circ\in F\}$. For integer endpoints it uses the fact that $n^\circ$ is an integer table of $F$, which is valid with exact integer margins.
* $F$ may be random, because identity-row right-hand sides and the thresholds $Y^{NM}$, $Y^{ER}$ depend on data. The inclusions are per replication, so this does no harm.
* (i) and (ii) for identity rows follow because $P(n^\circ\notin F)=0$. The belief-row readings are exactly the two inequalities, and the conditional-level caveat is stated and reproduced by X6.
* Remark 4.2 is now consistent: $\alpha_C$ plays the role of $P(n^\circ\notin F)$.
* C8's two printed texts match the theorem. The register row R27 states the restriction.

**One residual error, N1 (C2(b), one phrase; error in a parenthetical).**
* The corrected (b) says the integer gap set is "the lattice inside the continuous interval" *"without unit rows, or with integer boxes"*. Proof step 3 says the same in the form "only when … without fractional boxes on several units".
* Fractional **totals** among C2(b)'s own unit rows (class (i) totals on $(r,c)$) have the same effect as fractional boxes. RC11: four units, no boxes, totals $y_0+y_1\ge0.2$ and $y_2+y_3\ge0.2$, and a between sign $Y\ge0.8$. The continuous range $[0.8,4]$ contains 1, but the integer values are $\{2,3,4\}$.
* Fix: replace "or with integer boxes" by "or when every unit row (box, known cell, total) has an integer right-hand side". In proof step 3, replace "only when … i.e. without fractional boxes on several units" by "when every unit row has an integer right-hand side; otherwise it can be strictly smaller (T5, S5)". ("Only when" also overstates: a dominating threshold can restore the full lattice.)
* Integer endpoints, R19(b) and the formula $[\max(Y^-_{\rm rd},\lceil\cdot\rceil),\min(Y^+_{\rm rd},\lfloor\cdot\rfloor)]$ are unaffected; only the description of when the full lattice appears is wrong.

Nothing else new. I also rechecked these wording changes and found them correct:
* the R19(b) item 3 clause "(for (iii), in the reduced system)";
* the Corollary 2.2 proviso;
* the R19(c)(3) rewording;
* the §4.1 split of identity rows and belief rows;
* the C4 and C8 contract text;
* the §9 items.

### 7.3 Can the register rows move from "proof review pending" to "proof reviewed"?

* **R27: yes.** All issues assigned to it (2, 3, 4, 10, 11) are resolved, and the restated union bounds and Lemma 4.1(b) are correct. Its status remains inherited from the model interval (the §4.3 table). (iii-b) remains a computed example, not P.
* **R19: yes, once the one-phrase fix N1 is made in Prop. C2(b) and its proof step 3.** The fix is textual and needs no further proof review; the coordinator can confirm it by reading the two sentences. Every other R19 issue (1, 5's main sentence, 6, 7, 8, 9, 12–15) is resolved.

---

## 8. Round 3: the added §4.3 status rows (6a at `claude/mcount-theory` bbf8cf46)

*I merged bbf8cf46 into this branch as a fast-forward. It also carries 6e0a684, the round-2 fix N1, which I checked: C2(b) and its proof step 3 now say "when every unit row (box, known cell, total) has an integer right-hand side", exactly as asked. `verify_phase6a.py` is unchanged since ed709f11 (36/36). I checked the five new §4.3 rows and the R27 register row against 2b (Thm 2, D6, D7), P3 (§0.5 table, Thm P3.5(a), Prop. P3.5(b2)(ii), Thm P3.8) and PHASE5a.*

| New row | Verdict | Reasoning checked |
|---|---|---|
| L6, and L5 with one group (point-identified) | **correct** | $\tilde B_{rc}=\sum_iw_{ir}\alpha_{c,r}$ is one linear functional of the no-intercept Goodman regression, which has fixed dimension $R(C-1)$ (2b D7(a); with one group L5 is L6). With an empty constraint system, P3.5(a)'s event $\mathcal E$ consists of the two-sided constant term alone. So $\hat k=z_{1-\alpha_{\rm est}/2}$, and steps 4–7 add $S_\varphi$; P3.8 with $t_E(\varphi)$ works the same way. There are no districts, so the regime does not matter. The model is P3-S1a with $m$ constant. |
| L5 with several groups | **correct, and correctly restricted to R1** | $\dim\beta=D^\star R(C-1)$ grows with the number of groups. P3.5(a) step 1 needs a fixed-dimension asymptotically normal $\hat\beta$, which holds in R1 only. In R2 each group's regression rests on a bounded number of units, so neither consistency nor HC3 variance estimation is available. The row's "R2: simulation finding only" matches P3 §0.5 for asymptotic intervals. P3.8 still gives a finite-sample guarantee in R2 (wide); see nit (b). |
| pooled L4 | **correct, including the near-zero-width case** | 2b D7(b): $K_{c,rs}=Q_{c,rr}-Q_{c,rs}$ is affine in $\beta$ with fixed coefficients on $k$, so each equality becomes two rows $\pm(s(\beta)+Ak)\ge0$, as with R10's equalities. P3.5(a)'s proof needs only $k_0\in\hat\Theta^+$ on $\mathcal E$ (step 2). It does not need Slater or an interior point, so it holds when $\Theta$ is a single point (2×2, where L4 pins every coordinate). Since $\mathcal E$ contains both one-sided terms of each equality, the equality is in effect controlled two-sided. (b2)(ii) is one simultaneous event for all functionals. P3 §5's "width near zero" caveat applies to per-endpoint one-sided (CLR-type) intervals, not to these. The rung is pooled and fixed-dimensional, so it holds in any regime. |
| pooled M1(δ), one group | **correct** | The rows are $\lVert k^{(d)}-\bar k-c_d\rVert_\infty\le\delta$ with $\bar k$ free (2b §1 table; D6(d)). With one group, $c_1=0$ ($\zeta_{c,1}=0$) and $\bar k$ absorbs $k^{(1)}$, so the restriction is vacuous and the rung is M0. |
| grouped L4, M1(δ) without sign rows | **correct, R1 only** | The constraints stay of P3.5(a)'s form: $c_d$ are coordinates of $\operatorname{skew}(\zeta_{\cdot,d}\mathbf 1^\top)$ with $\zeta_{c,d}=\operatorname{diag}(Q^{(d)}_c-Q^{(1)}_c)$, affine in $\beta$. δ enters only as a constant. The dimension is fixed only with fixed groups, which matches P3's table for M1/M2 and 2c M5's scope sentence. Grouped rungs **with** sign rows remain "no theorem" (R16), and that row is unchanged. |

**R27 register row:** consistent with the table. It lists pooled M0, L4 and M1(δ) in any regime; M1/M2 and grouped L4/M1(δ) in R1; L6 and one-group L5 in any regime; multi-group L5 in R1; P3.8 finite-sample; CLR algorithmic; grouped rungs with signs none.

**Nits (wording; neither blocks):**
* (a) Some constraint terms do not depend on $\beta$, for example the δ-rows of the reference district, where $c_1=0$. Their $\hat\sigma=0$, and P3.5(a)'s max of standardised terms should be taken over the $\beta$-dependent terms only, with the rest held exactly. The proof is unaffected because these rows have no estimation error. (b2)(ii) handles this automatically. A one-line note under the table, or in the code contract, would prevent a division by zero.
* (b) The rows "L5 with several groups" and "grouped L4, M1(δ)" could add "R2: P3.8 (finite sample, wide)", so that "simulation finding only" is not read as "nothing proved in R2". The P3.8 row already covers it.

**Verdict: the five rows and the R27 register update are correct. The L4 near-zero-width case is covered, and multi-group L5 is correctly restricted to R1. R27 keeps "proof reviewed".**

*Addendum (6a at `claude/mcount-theory` 79e735b, merged here).* Both round-3 nits are fixed correctly.
* (a) A note under the §4.3 table now says that $\beta$-free constraint terms are held exactly, and that $\hat k$'s maximum is taken over the $\beta$-dependent terms and the constant only.
* (b) The multi-group L5, grouped L4/M1(δ) and grouped-R2 rows now cite P3.8 (finite sample, wide) for R2. That use is valid under P3.8's own hypotheses (P3-S1a, a correct mean model, a full-rank regression design).

Nothing else changed. **R19 and R27 remain "proof reviewed".**
