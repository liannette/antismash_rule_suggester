# Non-core segment groups: annotation review

2026-11-15. Segment grouping (step 1) on `data/fasta` (66 building blocks, 112 subclusters), with the
`absorb.py` fix applied. Of the 24 building blocks with more than one subcluster, 16 are all-core
and 8 have non-core groups. Those 8 are listed below.

For every subcluster missing a group, the group's proteins were searched with DIAMOND (blastp,
`--more-sensitive`, e < 1e-5) against **all** CDSs of that subcluster's MIBiG record:

- **no hit:** the record has no homolog, so the gene is either really absent or wrongly added to the other subcluster(s);
- **hit already in the gene list:** the homolog is curated but was placed in another group (low similarity or a paralog).

No search found a homolog that is missing from `genes_mibig`. All four cases flagged below were
checked against the records and the cited papers (see **Resolutions** at the end): all are real
biology, and the only gene-level change went the other way — `inkH` was removed from SC0070.

## Checked (see Resolutions below)

### BB0003 L-dihydrostreptose → L-rhamnose

| SC | BGC | product | G01 | G02 | G03 | G04 | G05 |
|---|---|---|---|---|---|---|---|
| SC0003 | BGC0000724 | streptomycin | strD | strE | strL | **strM** | — |
| SC0070 | BGC0000813 | K-252a | inkA | inkB | inkL | — | ~~inkH~~ |

**Resolved.** inkH was not the epimerase and has been removed from SC0070; the block was renamed
L-dihydrostreptose → **L-rhamnose**. K-252a has no epimerase at all. See Resolutions.

### BB0023 L-rhodinose / BB0025 D-olivose (both landomycin, BGC0000239)

| SC | BGC | product | G01 | G02 | G03 | G04 | G05 | G06 | G07 | G08 |
|---|---|---|---|---|---|---|---|---|---|---|
| SC0047 | BGC0001046 | streptolydigin | slgS1 | slgS2 | slgS3 | slgS4 | slgS6 | slgS5 | slgS7 | — |
| SC0148 | BGC0000239 | landomycin | lanG | lanH14 | lanS | lanT | lanQ | lanZ3 | lanZ1 | **lanZ2** |

| SC | BGC | product | G01 | G02 | G03 | G04 | G05 | G06 |
|---|---|---|---|---|---|---|---|---|
| SC0036 | BGC0000239 | landomycin | lanG | lanH14 | lanR | lanS | lanT | **lanZ2** |
| SC0052 | BGC0000036 | chlorothricin | chlC1 | chlC2 | chlC5 | chlC3 | chlC4 | — |

Both landomycin rows list two NDP-hexose synthetases (lanG and lanZ2). The only hit for lanZ2 in
the other record is its lanG counterpart (slgS1 41 % id; chlC1 43 % id), so the other clusters have
just one.

**Resolved.** lanZ2 belongs: landomycin genuinely has two synthetases, and lanZ2 is
landomycin-A-specific. See Resolutions.

### BB0007 3-amino-4,7-dihydroxycoumarin

| SC | BGC | product | G01 | G02 | G03 | G04 |
|---|---|---|---|---|---|---|
| SC0007 | BGC0000832 | clorobiocin | cloH | cloI | cloJ | cloK |
| SC0011 | BGC0000833 | coumermycin A1 | couH | couI | couJ | couK |
| SC0015 | BGC0000834 | novobiocin | novH | novI | novJ | novK |
| SC0039 | BGC0001072 | simocyclinone | simD6 | simD1 | simD2 | simD3 |
| SC0043 | BGC0000141 | rubradirin | rubC1 | rubC2 | rubC3 | — |

Rubradirin has no cloK/couK/novK/simD3 counterpart. Its only hit is rubC3 (36 % id to simD3), which
is already grouped with cloJ.

**Resolved.** Neither: NovJ/NovK are two subunits of one oxidase and the rub cluster has no NovK
homolog, as the source paper states. rubC1 is also a NovL+NovH domain fusion. See Resolutions.

### BB0018 3-amino-5-hydroxybenzoic acid (AHBA)

| SC | BGC | product | G01–G06 | G07 | G08 |
|---|---|---|---|---|---|
| SC0028 | BGC0000090 | macbecin | mbcO, mbcE, mbcD, mbcC, mbcL, mbcQ | — | — |
| SC0030 | BGC0000106 | naphthomycin | natG, natK, natL, natM, natJ, natN | natH | natI |
| SC0040 | BGC0000141 | rubradirin | rubG, rubK, rubL, rubM, rubJ, rubN | rubH | — |
| SC0044 | BGC0000136 | rifamycin | rifG, rifK, rifL, rifM, rifJ, rifN | rifH | rifI |

There is no rifH homolog in BGC0000090 and no rifI homolog in BGC0000090 or BGC0000141.

**Resolved.** Real biology, confirmed by a three-query search against all 23 macbecin CDSs. rifI is
auxiliary; the aDAHP synthase is supplied from outside the cluster. See Resolutions.

## Explained: real differences, not annotation errors

- **BB0059 4-methyl-3-hydroxyanthranilic acid:** ABW71854.1 (kynurenine 3-monooxygenase) is only in
  anthramycin (SC0169). The SC0095 actinomycin notes already say the cluster has no KMO.
- **BB0062 salicylic acid:** pyochelin uses pchA + pchB, while mycobactin uses the bifunctional
  salicylate synthase mbtI. pchA and mbtI are homologous (30 % id over about 55 %) but too divergent
  to group, and pchB has no counterpart in mycobactin. (Mechanism refined — see Resolutions.)
- **BB0088 2-deoxystreptamine:** btrE (butirosin, G04) and rbmC/neoE/parE/livE (G03) are the same
  dehydrogenase step, but only about 22 % identical, so they split into two groups. Every subcluster has one.

## Resolved by the absorb.py fix

These building blocks are all-core now. Before the fix they showed spurious non-core groups:
BB0016, BB0020, BB0073. BB0003, BB0007, BB0018, BB0023 and BB0025 also lost spurious groups.

## Resolutions (2026-11-15)

The four "worth checking" cases above were checked against the MIBiG records and the cited papers.
All four are real biology; no `genes_mibig` entry was wrong, and one gene was removed.

### BB0003 — resolved, and the building block was renamed

`inkH` is **not** the epimerase and has been dropped from SC0070 (`inkA;inkB;inkL`). Fig. 4 of
PMID 17396254, the paper's own pathway, contains only InkA, InkB, InkL and InkM; inkH appears
nowhere on it, and its Table 1 "Epimerase" label rests on 29 % id / 43 % sim to OvmX. A strM search
against all 14 CDSs of BGC0000813 gives no DIAMOND hit (ultra-sensitive, e ≤ 10) and local-alignment
z-scores of −1.1 to 1.2 against shuffled controls, so no gene in the cluster does this step.

The block itself was wrong: Table 1 calls InkL's homolog "dTDP-dihydrostreptose synthase StrL", but
StrL is a dTDP-4-dehydrorhamnose **reductase** (PMID 17390001, which calls it "StrL (reductase)";
MIBiG BGC0000724 annotates it "dTDP-4-dehydrorhamnose reductase"). PMID 17390001 also states that
all the enzymatic functions needed for dTDP-L-dihydrostreptose are predictable from the cluster
"except for the conversion of dTDP-L-rhamnose (33) to 34", leaving that step unassigned. Both subclusters therefore make **dTDP-L-rhamnose**, and BB0003 was renamed from
L-dihydrostreptose to L-rhamnose (PubChem CID 25310, InChIKey SHZGCJCMOBCMKK-JFNONXLTSA-N).

**Rule:** K-252a is a 3-gene cassette with no epimerase; streptomycin has the full 4. Make the
RmlC-type epimerase optional.

### BB0023 / BB0025 — lanZ2 belongs

Landomycin really has two NDP-hexose synthetases. PMID 22186970: "The genes lanG and lanZ2 encode
NDP-hexose synthases… These three enzymes [with lanH] catalyze the conversion of D-glucose-1-phosphate
to TDP-4-keto-6-deoxy-D-glucose." PMID 9933932 gives LanG 69 % / LanZ2 46 % id to GraD; MIBiG
annotates both "NDP-hexose synthetase homolog". lanZ2 is landomycin-A-specific — one of three ORFs
(with lanK and lanGT3) present in the *lan* cluster but absent from the *lnd* cluster of
*S. globisporus*, matching landomycin A's hexasaccharide versus landomycin E's trisaccharide.

**Rule:** allow 1–2 copies of the NDP-hexose synthetase step.

### BB0007 — rubradirin is missing a subunit, not a step

PMID 18080113 states it directly: "RubC3 is similar to NovJ, which, together with NovK, has been
loosely implicated in the conversion of hydroxytyrosine to keto-tyrosine, although **no NovK homologue
was detected in the rub gene cluster**". No NovC homolog either, and the authors note that some
4-hydroxyaminocoumarin genes "have yet to be discovered". NovJ and NovK are two subunits of one
oxidase (PMID 16171397; Pacholec et al. 2005, *Biochemistry* 44:12819-26), so this is a missing subunit.

`rubC1` is also a fusion of the NovL amide-synthetase and NovH tyrosine-adenylation domains in one
1317-aa ORF, so the cassette does not map 1:1 onto novHIJK. It is kept despite carrying the NovL
exit function (amide bond formation — excluded as novL/cloL/couL are elsewhere in BB0007) because
the same ORF performs the NovH loading step.

**Rule:** fourth gene optional; expect rubC1 to hit two domain profiles.

### BB0018 — macbecin genuinely lacks both genes

Confirmed by search, not just by annotation: rifH, natH and rubH give no DIAMOND hit against any of
the 23 CDSs of BGC0000090 (all local-alignment z ≤ 2.7). For rifI the best candidate mbcU scores 55
(41 % id over 71 aa, empirical p ≈ 0.013) against 942 for the true rifI–natI pair — best-of-noise.
mbcH/mbcI share only the gene letters: they are the methoxymalonyl-ACP cassette (SC0029/BB0020,
matching fkbG–fkbK), confirmed by Figure 2 of PMID 18800759, where they are coloured blue and the
AHBA genes purple.

rifI/natI are shikimate/quinate dehydrogenases, auxiliary to the committed route, which explains
their absence from rubradirin too. The aDAHP synthase is the first committed step, so its absence
from macbecin matters more; PMID 18800759 nonetheless calls the cluster "complete" for AHBA. Note
the producer *A. pretiosum* ATCC 31280 is the same strain as the ansamitocin producer (BGC0001511),
whose AHBA genes are themselves dispersed outside the MIBiG boundary (SC0200), so a shared genomic
AHBA supply is plausible but unverified. SC0028 was kept as accepted, consistent with SC0040,
SC0043 and SC0070, which also have documented missing steps.

**Rule:** six steps are present in all four subclusters and form the core requirement — aDHQ synthase
(rifG/natG/rubG/mbcO), aDHQ dehydratase (rifJ/natJ/rubJ/mbcL), AHBA synthase and oxidoreductase
(rifK+rifL and counterparts, mbcD+mbcE), phosphatase (rifM/natM/rubM/mbcC) and kinase
(rifN/natN/rubN/mbcQ). The two optional ones are the aDAHP synthase (3/4: rifH, natH, rubH; absent
from macbecin) and the shikimate/quinate dehydrogenase (2/4: rifI, natI only).

### BB0062 — correction to the note above

MbtI is bifunctional, but the mechanism is weaker than "bifunctional" suggests. PMID 17240979
(Zwahlen et al., *Biochemistry* 46:954, already SC0195's reference) shows isochorismate is a
kinetically competent intermediate — below pH 7.5 it is the dominant product, above it salicylate
forms without isochorismate accumulating. MbtI belongs to the MST family and **lacks** the active-site
residues unique to pyruvate-eliminating members of it; without Mg²⁺ it shows promiscuous chorismate
mutase activity resembling PchB's, so the authors *propose* a sigmatropic elimination like PchB's.
pchB has no sequence counterpart in mbtI, so the grouping is correct as it stands.

**Rule:** the salicylate cassette is 1 or 2 genes. Requiring both would miss mycobactin.
