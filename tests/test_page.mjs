/* Offline tests for the page logic around the variant resolver.

   tests/test_hgvs.mjs covers the resolver's arithmetic. This covers the glue
   that CI otherwise would not: geneEntryFor, which decides *which* gene a
   variant refers to and holds the transcript-safety rules, and the rendering
   in variantPanel/showVariant.

   The whole <script> from site/index.html runs against a stubbed DOM and a
   stubbed fetch serving synthetic genes, so there are no fixture files and no
   network. Run: node tests/test_page.mjs                                   */

import fs from 'fs';
import path from 'path';

/* ---- synthetic data ----------------------------------------------------
   Two 3-exon genes, 30 bp CDS, one per strand, mirroring test_hgvs.mjs.
   PLUSG has exons at 1000-1009 / 2000-2009 / 3000-3009 on chr7; MINUSG the
   same coordinates on chr8 with the ranks reversed. One amplicon each, over
   exon 2 with 50 bp flanks, sequence in transcript orientation.          */
const CDS = 'ATGAAACCCGGGTTTAAACCCGGGTTTTAA';

function mkGene(name, chrom, strand, transcript, refseq){
  const base = [
    {n:1, start:1000, end:1009, cds_start:1000, cds_end:1009, coding:true, len:10},
    {n:2, start:2000, end:2009, cds_start:2000, cds_end:2009, coding:true, len:10},
    {n:3, start:3000, end:3009, cds_start:3000, cds_end:3009, coding:true, len:10},
  ];
  const exons = strand === 1 ? base : base.map((e,i)=>({...base[base.length-1-i], n:e.n}));
  const thirds = [CDS.slice(0,10), CDS.slice(10,20), CDS.slice(20,30)];
  const start = 1950, end = 2059;
  let seq = '';
  for(let p = start; p <= end; p++){
    if(p < 2000 || p > 2009){ seq += 'n'; continue; }
    const ex2 = exons.find(e=>e.n===2);
    seq += thirds[1][strand===1 ? p-ex2.cds_start : ex2.cds_end-p];
  }
  if(strand === -1) seq = seq.split('').reverse().join('');
  return {
    gene:name, gene_id:'ENSG00000000001', transcript, refseq, selection:'MANE Select',
    assembly:'GRCh38', chrom, strand, exon_count:3,
    annotated:{snps:true,repeats:true,segdups:true}, exons, n_snps:0,
    amplicons:[{label:'2', exons:[2], fragment:'', target_start:2000, target_end:2009,
      target_len:10, relaxed:false, tier:'clean', status:'OK', segdups:[], warnings:[],
      pairs:[{rank:1, f_seq:'ACGTACGTACGTACGTAC', r_seq:'TGCATGCATGCATGCATG',
        f_tm:60, r_tm:60, f_gc:50, r_gc:50, f_len:18, r_len:18, product:end-start+1,
        penalty:1, f_dist:50, r_dist:50, chrom, start, end, amplicon:seq,
        f_snps:[], r_snps:[], f_repeat:0, r_repeat:0,
        f_tm_na:54.2, r_tm_na:54.0, d_tm:0,
        ispcr:{checked:true,on_target:true,off_target:[]}}]}],
    failures:[],
  };
}

const PLUSG  = mkGene('PLUSG',  '7',  1, 'ENST00000000001.3', 'NM_000001.4');
/* Rank 1 mismatched (dTm 4), rank 2 matched (dTm 0.5): lets the dTm
   preference be exercised without touching the other genes. */
{
  const a = PLUSG.amplicons[0], p1 = a.pairs[0];
  p1.f_tm = 62; p1.r_tm = 58;
  p1.d_tm = 4;
  a.pairs.push({...p1, rank:2, f_tm:60, r_tm:60.5, d_tm:0.5, f_seq:'ACGTACGTACGTACGTAA',
                r_seq:'TGCATGCATGCATGCATT', penalty:2});
}
const MINUSG = mkGene('MINUSG', '8', -1, 'ENST00000000002.1', 'NM_000002.2');

const INDEX = {
  assembly:'GRCh38', ucsc_db:'hg38', generated:'2026-01-01', version:'2.1.0',
  params:{min_dist:50, primer_region:150, product_min:200, product_max:500,
    max_target:320, overlap:50, combine:true, target:'coding', tm_min:57, tm_opt:60,
    tm_max:63, size_min:18, size_opt:20, size_max:27, gc_min:30, gc_max:70,
    gc_clamp:1, snp_maf:0.01, snp_3prime:5, num_return:3, relax:true},
  ispcr:true, tiers:{}, errors:[{gene:'BROKEN', error:'lookup failed'}],
  genes:[
    {gene:'PLUSG',  transcript:PLUSG.transcript,  refseq:PLUSG.refseq,  chrom:'7',
     amplicons:1, ok:1, check:0, failures:0, start:1000, end:3009},
    {gene:'MINUSG', transcript:MINUSG.transcript, refseq:MINUSG.refseq, chrom:'8',
     amplicons:1, ok:1, check:0, failures:0, start:1000, end:3009},
  ],
};
const FILES = {
  'data/index.json': INDEX,
  'data/genes/PLUSG.json': PLUSG,
  'data/genes/MINUSG.json': MINUSG,
};

/* ---- stub DOM + fetch, then run the real page script ---- */
const els = {};
const el = id => (els[id] = els[id] || {id, innerHTML:'', textContent:'', value:'',
  classList:{add(){},remove(){}}, addEventListener(){}, dataset:{}, onclick:null});
globalThis.document = {
  querySelector: s => el(s.replace('#','')),
  querySelectorAll: () => [],
};
globalThis.location = {hash:'', pathname:'/exonprimer-batch/', search:''};
globalThis.history = {pushState(){}};
globalThis.window = {addEventListener(){}};
Object.defineProperty(globalThis, 'navigator',
  {value:{clipboard:{writeText:()=>Promise.resolve()}}, configurable:true});
/* UCSC API stub: NM_777777.1 is a second transcript of PLUSG that numbers
   its CDS 9 bp later (an extra upstream coding exon), and shares exon 2
   exactly. UCSC rows are 0-based half-open, as the real API serves them. */
const UCSC_ROW = {
  name:'NM_777777.1', name2:'PLUSG', chrom:'chr7', strand:'+',
  txStart:899, txEnd:3009, cdsStart:899, cdsEnd:3009,
  exonStarts:'899,999,1999,2999,', exonEnds:'909,1009,2009,3009,',
};
let ucscCalls = 0;
globalThis.fetch = async u => {
  if(u in FILES) return {json: async () => JSON.parse(JSON.stringify(FILES[u]))};
  if(u.startsWith('https://api.genome.ucsc.edu/')){
    ucscCalls++;
    if(u.includes('/search')){
      if(!u.includes('NM_777777')) return {ok:true, json:async()=>({positionMatches:[]})};
      return {ok:true, json: async () => ({positionMatches:[{name:'ncbiRefSeqCurated',
        matches:[{position:'chr7:900-3009', posName:'PLUSG (NM_777777.1)'}]}]})};
    }
    if(u.includes('track=ncbiRefSeqCurated'))
      return {ok:true, json: async () => ({ncbiRefSeqCurated:[UCSC_ROW]})};
    return {ok:true, json: async () => ({})};
  }
  throw new Error('404 ' + u);
};

const page = fs.readFileSync(path.join(import.meta.dirname, '..', 'site', 'index.html'), 'utf8');
const script = page.match(/<script>([\s\S]*)<\/script>/);
if(!script){ console.error('no <script> block in site/index.html'); process.exit(1); }
const api = new Function(script[1] + `
  return {init, route, showVariant, geneEntryFor, parseHGVS, submitQuery,
          renderGene, geneTSV, ucscRowToModel, ucscTranscript, resolveViaUcsc,
          orderedPairs, bestPair, dTm, setPref:v=>{PREF_MATCHED=v},
          view:()=>document.querySelector('#view').innerHTML};`)();

let pass = 0, fail = 0;
const ok = (c, m) => { if(c) pass++; else { fail++; console.log('  FAIL: ' + m); } };
const view = () => api.view();
const has  = t => view().includes(t);

await api.init();

/* ---- init ---- */
ok(els['meta'].innerHTML.includes('2 genes'), 'header reports the gene count');
ok(els['genelist'].innerHTML.includes('PLUSG'), 'datalist is populated');
ok(els['note'].innerHTML.includes('Primer3'), 'method note is rendered');
ok(els['note'].innerHTML.includes('BROKEN'), 'genes that failed are listed');
ok(!els['note'].innerHTML.includes('undefined'), 'no undefined in the method note');

/* ---- geneEntryFor: which gene does this variant mean? ---- */
ok(api.geneEntryFor(api.parseHGVS('PLUSG:c.15T>G')).hit.gene === 'PLUSG', 'gene by symbol');
ok(api.geneEntryFor(api.parseHGVS('plusg:c.15T>G')).hit.gene === 'PLUSG', 'symbol is case-insensitive');
ok(api.geneEntryFor(api.parseHGVS('ENST00000000001:c.15T>G')).hit.gene === 'PLUSG', 'gene by Ensembl transcript');
ok(api.geneEntryFor(api.parseHGVS('NM_000001.4:c.15T>G')).hit.gene === 'PLUSG', 'gene by exact RefSeq');

const ver = api.geneEntryFor(api.parseHGVS('NM_000001.2:c.15T>G'));
ok(ver.hit && ver.hit.gene === 'PLUSG', 'older RefSeq version still resolves');
ok(ver.warn && /NM_000001\.4/.test(ver.warn), 'version difference is flagged, not silent');

const unk = api.geneEntryFor(api.parseHGVS('NM_999999.9:c.15T>G'));
ok(!unk.hit && unk.lookup === 'NM_999999.9',
   'unknown transcript is never renumbered against our own, it is offered for lookup');

const byPos = api.geneEntryFor(api.parseHGVS('chr7:2005'));
ok(byPos.hit && byPos.hit.gene === 'PLUSG', 'gene found from a coordinate');
ok(api.geneEntryFor(api.parseHGVS('chr7:999999')).error, 'coordinate outside every gene is refused');
ok(api.geneEntryFor(api.parseHGVS('chr22:2005')).error, 'coordinate on an unused chromosome is refused');

/* ---- rendering ---- */
await api.showVariant('PLUSG:c.15T>G');
ok(has('PLUSG_E2'), 'covering amplicon is listed');
ok(has('chr7:2,004') || has('2,004'), 'mapped coordinate is shown');
ok(has('reference base T confirmed'), 'reference base is confirmed');
ok(has('ACGTACGTACGTACGTAC'), 'forward primer is shown');
ok(has('TGCATGCATGCATGCATG'), 'reverse primer is shown');
ok(!has('undefined'), 'no undefined leaked into the panel');
ok(has('#PLUSG/2'), 'amplicon links back to the gene and exon');

await api.showVariant('MINUSG:c.15T>G');
ok(has('MINUSG_E2') && has('reference base T confirmed'), 'minus-strand gene resolves and confirms its base');

await api.showVariant('PLUSG:c.15A>G');
ok(has('Reference mismatch'), 'wrong reference allele is surfaced to the reader');

await api.showVariant('PLUSG:p.Phe5Leu');
ok(has('codon 5') && has('PLUSG_E2'), 'protein variant resolves to its codon and amplicon');

await api.showVariant('PLUSG:c.-45G>A');
ok(has('UTR'), 'UTR variant explains itself instead of showing an exon');

await api.showVariant('PLUSG:c.20+500G>A');
ok(has('No amplicon in this panel covers'), 'deep intronic says no amplicon covers it');
ok(has('150'), 'and says how far the primers actually reach');

await api.showVariant('PLUSG:c.99999A>G');
ok(has('past the end'), 'position beyond the CDS is refused');

await api.showVariant('utter nonsense');
ok(has('Could not read'), 'unparseable input gets a usable message');

/* ---- a real symbol outside the panel is not "unreadable HGVS" ---- */
await api.showVariant('TUBB4B');
ok(has('not in our database') && !has('Could not read'),
   'a symbol-shaped query absent from the panel says so instead of blaming HGVS');
ok(has('2 genes'), 'and says how many genes the panel covers');

await api.showVariant('PLUSH');
ok(has('Did you mean') && has('>PLUSG<'), 'a one-letter typo suggests the real symbol');

await api.showVariant('NM_999999.9');
ok(has('not in our database'), 'a bare accession is treated as a missing reference too');

await api.showVariant('TUBB4B:c.100G>A');
ok(has('not a gene, transcript or RefSeq in our database'),
   'and the same wording covers an HGVS variant on a gene we do not carry');

await api.showVariant('NM_999999.9:c.15T>G');
ok(has('is not the transcript this panel uses') && has('Resolve NM_999999.9 via UCSC'),
   'an unknown transcript offers the UCSC conversion instead of a dead end');
ok(!has('<table'), 'and renders no primers until the conversion is actually run');

await api.showVariant('chr7:2005');
ok(has('PLUSG_E2'), 'bare coordinate finds the amplicon');

/* ---- dTm: shown, and optionally preferred (A2/A3) ---- */
globalThis.location.hash = '';
api.setPref(false);
await api.showVariant('PLUSG:c.15T>G');
ok(has('\u0394Tm'), 'the variant table carries a dTm column');

const amp = JSON.parse(JSON.stringify(PLUSG)).amplicons[0];
ok(api.dTm(amp.pairs[0]).toFixed(1) === '4.0', 'dTm is read off the published Tm values');
ok(api.bestPair(amp) === amp.pairs[0], 'with the preference off, Primer3 rank 1 is shown');
api.setPref(true);
ok(api.bestPair(amp).rank === 2, 'with it on, the matched alternative comes first');
ok(api.orderedPairs(amp).length === amp.pairs.length, 'no pair is dropped by the re-ordering');
api.setPref(false);
ok(api.bestPair(amp).rank === 1, 'and the preference is reversible');

/* A promotion must never trade away genome-wide specificity. */
const risky = JSON.parse(JSON.stringify(PLUSG)).amplicons[0];
risky.pairs[1].ispcr = {checked:true, on_target:true, off_target:[{chrom:'chr9'}]};
api.setPref(true);
ok(api.bestPair(risky).rank === 1,
   'a matched alternative with an off-target product is not promoted');
ok(api.orderedPairs(risky).length === risky.pairs.length, 'and it is still listed');
api.setPref(false);

/* A matched rank 1 must never be displaced. */
const good = JSON.parse(JSON.stringify(PLUSG)).amplicons[0];
good.pairs[0].f_tm = 60; good.pairs[0].r_tm = 60; good.pairs[0].d_tm = 0;
api.setPref(true);
ok(api.bestPair(good).rank === 1, 'an already-matched rank 1 is left alone');
api.setPref(false);

/* The TSV export carries the delta too. */
ok(api.geneTSV(PLUSG).split('\n')[0].includes('dTm_pair'), 'TSV export includes dTm');

/* The published d_tm wins over recomputing it, so the three numbers a reader
   sees always agree; a build without it still works. */
const withD = JSON.parse(JSON.stringify(PLUSG)).amplicons[0];
ok(api.dTm(withD.pairs[0]) === 4, 'the published d_tm is used when present');
delete withD.pairs[0].d_tm;
ok(api.dTm(withD.pairs[0]) === 4, 'and is recomputed for a build that predates it');

await api.showVariant('PLUSG:c.15T>G');
ok(has('ΔTm'), 'the variant table still carries ΔTm with the published field');

/* ---- an unknown transcript resolved through UCSC (A4/A5) ---- */
const before = ucscCalls;
const txm = api.ucscRowToModel(UCSC_ROW, 'chr7');
ok(txm.gene === 'PLUSG' && txm.strand === 1, 'a UCSC row yields the symbol and strand');
ok(txm.exons.length === 4 && txm.exons[0].start === 900 && txm.exons[0].end === 909,
   'UCSC 0-based half-open blocks become our 1-based inclusive exons');
ok(txm.exons[2].start === 2000 && txm.exons[2].end === 2009,
   'and exon 3 of that transcript is exon 2 of ours, same coordinates');
ok(ucscCalls === before, 'building a model from a row costs no request');

/* c.25 on NM_777777.1 is c.15 on ours: one extra 10-bp coding exon upstream. */
const vTx = api.parseHGVS('NM_777777.1:c.25T>G');
const res = await api.resolveViaUcsc(vTx);
ok(!res.error && !res.notInPanel, 'the foreign transcript resolves');
ok(res.pos === 2004, `it maps to the same genomic base as PLUSG c.15 (got ${res.pos})`);
ok(res.gene.gene === 'PLUSG', 'and lands on our own gene record');
ok(res.amplicons.length === 1, 'the covering amplicon is found');
ok(res.notes.some(n=>/NM_777777\.1/.test(n)), 'the panel says which transcript did the numbering');
ok(res.notes.some(n=>/exon 3 of NM_777777\.1 is exon 2 of/.test(n)
                  && /same coordinates, so the same amplicon/.test(n)),
   'renumbered exon is matched by coordinates, so the same amplicon is offered');
ok(res.exon && res.exon.n === 2, 'and the exon shown is ours, not the foreign rank');
ok(res.warnings.some(w=>/numbering came from/.test(w)),
   'a warning keeps the numbering difference visible');
ok(res.notes.some(n=>/reference base T confirmed/.test(n)),
   'the reference base is checked against our own amplicon sequence');

/* A transcript of a gene outside the panel must not be renumbered. */
const away = api.parseHGVS('NM_888888.1:c.25T>G');
await api.showVariant('NM_888888.1:c.25T>G');
ok(has('Resolve NM_888888.1 via UCSC'), 'an unknown accession offers the lookup');
let threw = false;
try{ await api.resolveViaUcsc(away); }catch(e){ threw = /UCSC does not know/.test(e.message); }
ok(threw, 'and an accession UCSC cannot place raises a readable error');

/* ---- routing: a gene symbol must never be read as a variant ---- */
globalThis.location.hash = '';
api.submitQuery('PLUSG');
ok(location.hash === 'PLUSG', `gene symbol routes as a gene (got "${location.hash}")`);
api.submitQuery('PLUSG:c.15T>G');
ok(location.hash.startsWith('v='), `variant routes as a variant (got "${location.hash}")`);
api.submitQuery('  plusg  ');
ok(location.hash === 'PLUSG', 'whitespace and case are handled');

/* ---- a variant hash renders through route() ---- */
globalThis.location.hash = '#v=' + encodeURIComponent('PLUSG:c.15T>G');
await api.route();
ok(has('PLUSG_E2'), 'route() handles a v= hash');

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
