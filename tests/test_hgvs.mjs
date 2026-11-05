/* Offline tests for the HGVS lookup shipped in site/index.html.

   The resolver is extracted from the page itself, between the
   HGVS-RESOLVER markers, so these tests exercise the code that is served
   rather than a copy of it. Fixtures are synthetic genes on both strands,
   in the spirit of tests/test_offline.py: no network, no real data.

   Run: node tests/test_hgvs.mjs                                          */

import fs from 'fs';
import path from 'path';

const page = fs.readFileSync(path.join(import.meta.dirname, '..', 'site', 'index.html'), 'utf8');
const m = page.match(/\/\* HGVS-RESOLVER-START \*\/([\s\S]*?)\/\* HGVS-RESOLVER-END \*\//);
if (!m) { console.error('resolver markers not found in site/index.html'); process.exit(1); }
const R = new Function(m[1] + `
  return {parseHGVS, cdsToGenomic, cdsLength, refBaseAt, resolveInGene, amplsCovering, codingExons};`)();

let pass = 0, fail = 0;
const ok = (c, msg) => { if (c) pass++; else { fail++; console.log('  FAIL: ' + msg); } };

/* ---- fixtures -----------------------------------------------------------
   Two 3-exon genes with a 30 bp CDS, one per strand. PLUS has exons at
   1000-1009, 2000-2009, 3000-3009; MINUS mirrors it so that c.1 is the
   highest coordinate. Each carries one amplicon spanning exon 2 plus 50 bp
   of flank, with a sequence in transcript orientation.                    */
const CDS = 'ATGAAACCCGGGTTTAAACCCGGGTTTTAA';            // 30 bp, ATG..TAA
const comp = s => s.replace(/[ACGT]/g, c => ({A:'T',C:'G',G:'C',T:'A'}[c]));
const rc = s => comp(s).split('').reverse().join('');

function mkGene(strand){
  const ex = [
    {n:1, start:1000, end:1009, cds_start:1000, cds_end:1009, coding:true, len:10},
    {n:2, start:2000, end:2009, cds_start:2000, cds_end:2009, coding:true, len:10},
    {n:3, start:3000, end:3009, cds_start:3000, cds_end:3009, coding:true, len:10},
  ];
  // On the minus strand, exon rank 1 must be the highest coordinate.
  const exons = strand === 1 ? ex : ex.map((e,i) => ({...ex[ex.length-1-i], n:e.n}));
  // CDS thirds, per exon in transcript order
  const thirds = [CDS.slice(0,10), CDS.slice(10,20), CDS.slice(20,30)];
  // Amplicon over exon 2 (genomic 2000-2009) with 50 bp flanks.
  const start = 1950, end = 2059;
  let seq = '';
  for (let p = start; p <= end; p++){
    const inExon = p >= 2000 && p <= 2009;
    if (!inExon) { seq += 'n'; continue; }
    const ex2 = exons.find(e => e.n === 2);
    const off = strand === 1 ? p - ex2.cds_start : ex2.cds_end - p;
    seq += thirds[1][off];                       // transcript-orientation base
  }
  if (strand === -1) seq = seq.split('').reverse().join('');   // string runs end->start
  return {
    gene:'TESTG', transcript:'ENST00000000001.1', refseq:'NM_000001.1',
    chrom:'7', strand, selection:'MANE Select', exons,
    amplicons:[{label:'2', exons:[2], status:'OK', target_start:2000, target_end:2009,
      pairs:[{f_seq:'AAAA', r_seq:'TTTT', product:end-start+1, start, end, amplicon:seq}]}],
    failures:[],
  };
}
const PLUS = mkGene(1), MINUS = mkGene(-1);

/* ---- parser ---- */
const cases = [
  ['NM_000110.4:c.1129G>A', v => v.kind==='c' && v.cpos===1129 && v.refBase==='G' && v.altBase==='A' && v.ref==='NM_000110.4'],
  ['POLG:c.2243G>C',        v => v.cpos===2243 && v.ref==='POLG'],
  ['c.1129+5G>A',           v => v.cpos===1129 && v.offset===5],
  ['c.2243-12T>C',          v => v.offset===-12],
  ['c.-45G>A',              v => v.utr==='-'],
  ['c.*120A>G',             v => v.utr==='*'],
  ['p.Ala467Thr',           v => v.kind==='p' && v.codon===467 && v.aa==='A'],
  ['p.Trp748Ser',           v => v.codon===748 && v.aa==='W'],
  ['p.W748S',               v => v.codon===748 && v.aa==='W'],
  ['chr15:89327201C>T',     v => v.kind==='g' && v.chrom==='15' && v.pos===89327201],
  ['15:89327201',           v => v.kind==='g' && v.pos===89327201],
  ['NC_000015.10:g.89327201C>T', v => v.kind==='g' && v.pos===89327201],
  ['c.1129_1130del',        v => v.cpos===1129 && v.cpos2===1130],
  ['c.1129dup',             v => v.cpos===1129],
];
for (const [txt, check] of cases){
  const v = R.parseHGVS(txt);
  ok(v && check(v), `parse ${txt}`);
}
ok(R.parseHGVS('not a variant') === null, 'reject free text');
ok(R.parseHGVS('') === null, 'reject empty');

/* ---- CDS mapping, both strands ---- */
for (const [name, g] of [['plus', PLUS], ['minus', MINUS]]){
  ok(R.cdsLength(g) === 30, `${name}: CDS length 30`);
  const expect = g.strand === 1
    ? {1:1000, 10:1009, 11:2000, 20:2009, 21:3000, 30:3009}
    : {1:3009, 10:3000, 11:2009, 20:2000, 21:1009, 30:1000};
  for (const [c, want] of Object.entries(expect)){
    const r = R.cdsToGenomic(g, +c);
    ok(r && r.pos === want, `${name}: c.${c} -> ${want} (got ${r && r.pos})`);
  }
  ok(R.cdsToGenomic(g, 31) === null, `${name}: c.31 out of range`);
}

/* ---- reference base read-back, both strands ---- */
for (const [name, g] of [['plus', PLUS], ['minus', MINUS]]){
  let good = 0;
  for (let c = 11; c <= 20; c++){                       // exon 2, the amplicon
    const r = R.cdsToGenomic(g, c);
    if (R.refBaseAt(g, r.pos) === CDS[c-1]) good++;
  }
  ok(good === 10, `${name}: all 10 exon-2 reference bases read correctly (${good}/10)`);
}

/* ---- intron offsets step the right way per strand ---- */
{
  const p = R.resolveInGene(PLUS, R.parseHGVS('c.20+5G>A'));
  ok(!p.error && p.pos === 2014, `plus: c.20+5 -> 2014 (got ${p.pos})`);
  const m = R.resolveInGene(MINUS, R.parseHGVS('c.20+5G>A'));
  ok(!m.error && m.pos === 1995, `minus: c.20+5 -> 1995 (got ${m.pos})`);
  const p2 = R.resolveInGene(PLUS, R.parseHGVS('c.11-5G>A'));
  ok(!p2.error && p2.pos === 1995, `plus: c.11-5 -> 1995 (got ${p2.pos})`);
}

/* ---- resolution, refusals and the reference check ---- */
{
  const r = R.resolveInGene(PLUS, R.parseHGVS('c.15T>G'));    // c.15 is T
  ok(!r.error && r.amplicons.length === 1, 'exonic variant finds its amplicon');
  ok(r.notes.some(n => /reference base/.test(n)), 'correct reference allele is confirmed');

  const bad = R.resolveInGene(PLUS, R.parseHGVS('c.15A>G'));   // real base is T
  ok(bad.warnings.some(w => /Reference mismatch/.test(w)), 'wrong reference allele is flagged');

  ok(/UTR/.test(R.resolveInGene(PLUS, R.parseHGVS('c.-45G>A')).error || ''), "5'UTR refused, not mapped to an exon");
  ok(/UTR/.test(R.resolveInGene(PLUS, R.parseHGVS('c.*99A>G')).error || ''), "3'UTR refused");
  ok(/past the end/.test(R.resolveInGene(PLUS, R.parseHGVS('c.9999A>G')).error || ''), 'beyond CDS refused');
  ok(/chr9/.test(R.resolveInGene(PLUS, R.parseHGVS('chr9:1000A>G')).error || ''), 'wrong chromosome refused');

  const deep = R.resolveInGene(PLUS, R.parseHGVS('c.20+500G>A'));
  ok(!deep.error && deep.amplicons.length === 0, 'deep intronic maps but covers no amplicon');

  const prot = R.resolveInGene(PLUS, R.parseHGVS('p.Phe5Leu'));   // codon 5 = c.13-15 = TTT
  ok(!prot.error && prot.pos === 2002 && prot.posEnd === 2004, `p. maps to its codon (got ${prot.pos}-${prot.posEnd})`);

  const rng = R.resolveInGene(PLUS, R.parseHGVS('c.13_15del'));
  ok(!rng.error && rng.span && rng.span[0] === 2002 && rng.span[1] === 2004, 'range maps both endpoints');
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
