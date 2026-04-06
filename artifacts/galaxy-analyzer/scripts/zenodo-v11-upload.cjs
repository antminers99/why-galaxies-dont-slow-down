#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const LATEST_ID = '19434177';
const BASE = 'zenodo.org';

function zenodoRequest(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: urlPath,
      method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'User-Agent': 'GalaxyAnalyzer/11.0', 'Accept': 'application/json' },
    };
    if (body && contentType) {
      opts.headers['Content-Type'] = contentType;
      if (typeof body === 'string') opts.headers['Content-Length'] = Buffer.byteLength(body);
    }
    const req = https.request(opts, res => {
      let data = '';
      res.on('data', d => data += d);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('HTTP ' + res.statusCode + ': ' + data.substring(0, 500));
          reject(new Error('HTTP ' + res.statusCode));
        } else {
          try { resolve(JSON.parse(data)); } catch { resolve(data); }
        }
      });
    });
    req.on('error', reject);
    if (body && typeof body !== 'string') body.pipe(req);
    else { if (body) req.write(body); req.end(); }
  });
}

async function uploadFile(draftId, filePath, fileName) {
  const stat = fs.statSync(filePath);
  const initBody = JSON.stringify([{ key: fileName }]);
  await zenodoRequest('POST', '/api/records/' + draftId + '/draft/files', initBody, 'application/json');

  await new Promise((resolve, reject) => {
    const encodedName = encodeURIComponent(fileName);
    const opts = {
      hostname: BASE,
      path: '/api/records/' + draftId + '/draft/files/' + encodedName + '/content',
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'User-Agent': 'GalaxyAnalyzer/11.0',
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
      },
    };
    const req = https.request(opts, res => {
      let data = '';
      res.on('data', d => data += d);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('Upload error ' + res.statusCode + ' for ' + fileName + ': ' + data.substring(0, 300));
          reject(new Error('Upload content failed'));
        } else {
          resolve();
        }
      });
    });
    req.on('error', reject);
    fs.createReadStream(filePath).pipe(req);
  });

  await zenodoRequest('POST', '/api/records/' + draftId + '/draft/files/' + encodeURIComponent(fileName) + '/commit', '', 'application/json');
  console.log('  Uploaded: ' + fileName + ' (' + stat.size + ' bytes)');
}

async function main() {
  console.log('=== Zenodo v11 Upload: Beyond a0 — Decoding the Coupling Law ===');
  console.log('=== Phases 300-303: Physical Interpretation Program ===\n');

  let draftId;

  console.log('Step 1: Create new version from record ' + LATEST_ID + '...');
  try {
    const newVer = await zenodoRequest('POST', '/api/records/' + LATEST_ID + '/versions', '', 'application/json');
    draftId = newVer.id;
    console.log('  Created new draft ID: ' + draftId);
  } catch (e) {
    console.log('  Create failed (draft may already exist), checking latest draft...');
    const latestRec = await zenodoRequest('GET', '/api/records/' + LATEST_ID, null, null);
    if (latestRec.links && latestRec.links.latest_draft) {
      const draftUrl = latestRec.links.latest_draft;
      const draftInfo = await zenodoRequest('GET', new URL(draftUrl).pathname, null, null);
      draftId = draftInfo.id;
      console.log('  Found existing draft ID: ' + draftId);
    } else {
      throw new Error('Cannot find or create draft');
    }
  }

  console.log('\nStep 2: Check and delete existing files from draft...');
  try {
    const filesInfo = await zenodoRequest('GET', '/api/records/' + draftId + '/draft/files', null, null);
    const entries = filesInfo.entries || filesInfo;
    if (Array.isArray(entries) && entries.length > 0) {
      for (const f of entries) {
        try {
          await zenodoRequest('DELETE', '/api/records/' + draftId + '/draft/files/' + encodeURIComponent(f.key), null, null);
          console.log('  Deleted: ' + f.key);
        } catch (e) {
          console.log('  Could not delete ' + f.key + ': ' + e.message);
        }
      }
    } else {
      console.log('  No existing files to delete.');
    }
  } catch (e) {
    console.log('  Could not list files: ' + e.message);
  }

  console.log('\nStep 3: Upload new files...');
  const baseDir = path.join(__dirname, '..');
  const filesToUpload = [
    { local: 'scripts/phase126-m4-candidate.cjs', remote: 'scripts/phase126-m4-candidate.cjs' },
    { local: 'scripts/phase127-vflat-challenge.cjs', remote: 'scripts/phase127-vflat-challenge.cjs' },
    { local: 'scripts/phase128-decode-4th-axis.cjs', remote: 'scripts/phase128-decode-4th-axis.cjs' },
    { local: 'scripts/phase129-vflat-decomposition.cjs', remote: 'scripts/phase129-vflat-decomposition.cjs' },
    { local: 'scripts/phase130-split-4th-sector.cjs', remote: 'scripts/phase130-split-4th-sector.cjs' },
    { local: 'scripts/phase131-decode-vfresid.cjs', remote: 'scripts/phase131-decode-vfresid.cjs' },
    { local: 'scripts/phase132-vfresid-reducibility.cjs', remote: 'scripts/phase132-vfresid-reducibility.cjs' },
    { local: 'scripts/phase132a-halo-death-match.cjs', remote: 'scripts/phase132a-halo-death-match.cjs' },
    { local: 'scripts/phase132b-mediation-causal.cjs', remote: 'scripts/phase132b-mediation-causal.cjs' },
    { local: 'scripts/phase132c-external-robustness.cjs', remote: 'scripts/phase132c-external-robustness.cjs' },
    { local: 'scripts/phase133a-regime-law.cjs', remote: 'scripts/phase133a-regime-law.cjs' },
    { local: 'scripts/phase133b-second-channel.cjs', remote: 'scripts/phase133b-second-channel.cjs' },
    { local: 'scripts/phase133c-coupling-drivers.cjs', remote: 'scripts/phase133c-coupling-drivers.cjs' },
    { local: 'scripts/phase134-external-validation.cjs', remote: 'scripts/phase134-external-validation.cjs' },
    { local: 'external-validation/phase200-data-assembly.cjs', remote: 'scripts/phase200-data-assembly.cjs' },
    { local: 'external-validation/phase201-blind-prediction.cjs', remote: 'scripts/phase201-blind-prediction.cjs' },
    { local: 'external-validation/phase202-hierarchy-replication.cjs', remote: 'scripts/phase202-hierarchy-replication.cjs' },
    { local: 'external-validation/phase203-logMhost-improvement.cjs', remote: 'scripts/phase203-logMhost-improvement.cjs' },
    { local: 'external-validation/phase204-final-external-synthesis.cjs', remote: 'scripts/phase204-final-external-synthesis.cjs' },
    { local: 'external-validation/phase300-sample-salvage.cjs', remote: 'scripts/phase300-sample-salvage.cjs' },
    { local: 'external-validation/phase301-vfresid-drivers.cjs', remote: 'scripts/phase301-vfresid-drivers.cjs' },
    { local: 'external-validation/phase302-regime-law.cjs', remote: 'scripts/phase302-regime-law.cjs' },
    { local: 'external-validation/phase303-physical-interpretation.cjs', remote: 'scripts/phase303-physical-interpretation.cjs' },
    { local: 'external-validation/PROGRAM.md', remote: 'external-validation/PROGRAM.md' },
    { local: 'public/phase126-m4-candidate.json', remote: 'results/phase126-m4-candidate.json' },
    { local: 'public/phase127-vflat-challenge.json', remote: 'results/phase127-vflat-challenge.json' },
    { local: 'public/phase128-decode-4th-axis.json', remote: 'results/phase128-decode-4th-axis.json' },
    { local: 'public/phase129-vflat-decomposition.json', remote: 'results/phase129-vflat-decomposition.json' },
    { local: 'public/phase130-split-4th-sector.json', remote: 'results/phase130-split-4th-sector.json' },
    { local: 'public/phase131-decode-vfresid.json', remote: 'results/phase131-decode-vfresid.json' },
    { local: 'public/phase132-vfresid-reducibility.json', remote: 'results/phase132-vfresid-reducibility.json' },
    { local: 'public/phase132a-halo-death-match.json', remote: 'results/phase132a-halo-death-match.json' },
    { local: 'public/phase132b-mediation-causal.json', remote: 'results/phase132b-mediation-causal.json' },
    { local: 'public/phase132c-external-robustness.json', remote: 'results/phase132c-external-robustness.json' },
    { local: 'public/phase133a-regime-law.json', remote: 'results/phase133a-regime-law.json' },
    { local: 'public/phase133b-second-channel.json', remote: 'results/phase133b-second-channel.json' },
    { local: 'public/phase133c-coupling-drivers.json', remote: 'results/phase133c-coupling-drivers.json' },
    { local: 'public/phase134-external-validation.json', remote: 'results/phase134-external-validation.json' },
    { local: 'public/phase200-external-dataset.json', remote: 'results/phase200-external-dataset.json' },
    { local: 'public/phase201-blind-prediction.json', remote: 'results/phase201-blind-prediction.json' },
    { local: 'public/phase202-hierarchy-replication.json', remote: 'results/phase202-hierarchy-replication.json' },
    { local: 'public/phase203-logMhost-improvement.json', remote: 'results/phase203-logMhost-improvement.json' },
    { local: 'public/phase204-final-external-synthesis.json', remote: 'results/phase204-final-external-synthesis.json' },
    { local: 'public/phase300-sample-salvage.json', remote: 'results/phase300-sample-salvage.json' },
    { local: 'public/phase301-vfresid-drivers.json', remote: 'results/phase301-vfresid-drivers.json' },
    { local: 'public/phase302-regime-law.json', remote: 'results/phase302-regime-law.json' },
    { local: 'public/phase303-physical-interpretation.json', remote: 'results/phase303-physical-interpretation.json' },
    { local: 'public/stage-A-master-table.json', remote: 'data/stage-A-master-table.json' },
    { local: 'public/sparc-table.json', remote: 'data/sparc-table.json' },
    { local: 'public/phase58a2-tidal-expansion.json', remote: 'data/phase58a2-tidal-expansion.json' },
    { local: 'public/sparc-results.json', remote: 'data/sparc-results.json' },
    { local: 'public/replication/REPRODUCIBLE_RESULT.md', remote: 'REPRODUCIBLE_RESULT.md' },
    { local: 'public/replication/MANUSCRIPT-TRACK2.md', remote: 'MANUSCRIPT-TRACK2.md' },
  ];

  for (const f of filesToUpload) {
    const fullPath = path.join(baseDir, f.local);
    if (fs.existsSync(fullPath)) {
      await uploadFile(draftId, fullPath, f.remote);
    } else {
      console.log('  SKIPPED (not found): ' + f.local);
    }
  }

  console.log('\nStep 4: Update metadata...');
  const metadata = {
    metadata: {
      title: 'Structured Regime-Dependent Coupling Law for Per-Galaxy a0 Variation in SPARC Galaxies — with Physical Interpretation (v11 — Phases 126-134, 200-204, 300-303)',
      description:
        '<p>v11 — Complete scientific narrative: empirical discovery (Phases 126-134), external validation (Phases 200-204), and physical interpretation (Phases 300-303).</p>' +
        '<h3>Central Result</h3>' +
        '<p>Per-galaxy a<sub>0</sub> variation in SPARC is governed by a structured, regime-dependent baryon-halo coupling law. ' +
        'The dominant transferable channel (VfResid) encodes a <strong>regime-dependent dynamical integration effect</strong> — ' +
        'accumulated baryon-halo processing that genuinely amplifies in deeper potential wells. ' +
        'This is supported by external validation on N=94 independent galaxies and physical interpretation via hypothesis testing.</p>' +
        '<h3>What is new in v11 (Phases 300-303)</h3>' +
        '<ol>' +
        '<li><strong>Phase 300</strong>: Sample salvage recovers 35 additional external galaxies (N=59&rarr;N=94), preserving all correlations</li>' +
        '<li><strong>Phase 301</strong>: VfResid driver analysis — haloK is shared #1 predictor (r=0.598 int, 0.511 ext); pure structural models fail (LOO R&sup2;=&minus;0.144); irreducible fraction 35.8%/68.5% still predicts a<sub>0</sub></li>' +
        '<li><strong>Phase 302</strong>: Regime law — hidden residual (after removing haloK) predicts a<sub>0</sub> at r=0.299 low-V vs r=0.749 high-V internally; both slope and correlation strengthen (slope ratio 2.6:1), ruling out pure observability</li>' +
        '<li><strong>Phase 303</strong>: Physical interpretation — H4 (dynamical integration) wins (score 10) over H1 (halo response, 7), H3 (feedback, 7), H2 (assembly, 4); logMeanRun is the only regime-strengthening proxy; gas fraction fades at high-V</li>' +
        '</ol>' +
        '<h3>Established Results (Phases 126-204)</h3>' +
        '<ul>' +
        '<li>Hierarchical coupling law: Core (44.1%) + VfResid (+17pp) + lhOuter (+4pp) = 65.4% LOO gap</li>' +
        '<li>External validation: complete hierarchy replicates in N=59 independent galaxies, 8/8 checks pass, VfResid dominates by 30-76pp even after improving environmental proxy</li>' +
        '</ul>' +
        '<h3>Leading Physical Interpretation</h3>' +
        '<p>H4: Dynamical integration — an accumulated product of baryon-halo interaction that deepens with potential well depth. ' +
        'haloK captures the halo-amplitude component; logMeanRun provides a partial handle on integration depth. ' +
        'Gas fraction and morphology are relevant at low V<sub>flat</sub> but fade in the high-V<sub>flat</sub> regime where the coupling law is strongest, ' +
        'consistent with a transition from feedback-dominated to integration-dominated physics.</p>' +
        '<h3>Caveats</h3>' +
        '<ul>' +
        '<li>All data from SPARC survey; cross-survey replication needed</li>' +
        '<li>Hypothesis scoring is heuristic, not formal Bayesian model comparison</li>' +
        '<li>Hypotheses are not mutually exclusive</li>' +
        '<li>High-V<sub>flat</sub> subsamples remain small (N=8-21)</li>' +
        '<li>Not peer-reviewed</li>' +
        '</ul>' +
        '<p>Claim: structured, regime-dependent coupling law with external support and a leading physical interpretation based on dynamical integration.</p>' +
        '<p>All analysis scripts, external validation program, physical interpretation program, and machine-readable results included for full reproducibility. ' +
        'Data: SPARC (Lelli, McGaugh &amp; Schombert, 2016, AJ 152, 157).</p>',
      resource_type: { id: 'dataset' },
      publication_date: new Date().toISOString().split('T')[0],
      publisher: 'Zenodo',
      creators: [
        { person_or_org: { type: 'personal', given_name: 'Galaxy Rotation', family_name: 'Curve Analyzer' } }
      ],
      rights: [
        { id: 'cc-by-4.0' }
      ],
      subjects: [
        { subject: 'galaxy rotation curves' },
        { subject: 'SPARC' },
        { subject: 'radial acceleration relation' },
        { subject: 'dark matter' },
        { subject: 'baryon-halo coupling' },
        { subject: 'rotation curve analysis' },
        { subject: 'galaxy dynamics' },
        { subject: 'modified gravity' },
        { subject: 'MOND' },
        { subject: 'per-galaxy a0' },
        { subject: 'kinematic residual' },
        { subject: 'cross-validation' },
        { subject: 'external validation' },
        { subject: 'regime-dependent' },
        { subject: 'hierarchical coupling' },
        { subject: 'dynamical integration' },
        { subject: 'physical interpretation' },
        { subject: 'hypothesis testing' },
      ],
      additional_descriptions: [
        {
          description: 'v11 — Phases 126-134 + 200-204 + 300-303 (Complete Scientific Narrative with Physical Interpretation). ' +
            'Concept DOI: 10.5281/zenodo.19430633. ' +
            'Previous version: v10 (DOI 10.5281/zenodo.19434177). ' +
            'New in v11: P300 Sample Salvage (N=59 to N=94 external), ' +
            'P301 VfResid Drivers (haloK shared #1, pure structure fails, 35.8%/68.5% irreducible still predicts a0), ' +
            'P302 Regime Law (hidden residual activates with Vflat: r=0.299 low-V to 0.749 high-V, slope ratio 2.6:1), ' +
            'P303 Physical Interpretation (H4 dynamical integration score=10, H1/H3 tied at 7, H2 at 4; logFgas bug fixed). ' +
            'Updated MANUSCRIPT-TRACK2.md with Sections 3.9-3.12, revised Abstract/Section 4/5/6. ' +
            'Claim: structured, regime-dependent coupling law with external support and leading H4 dynamical integration interpretation.',
          type: { id: 'notes' }
        }
      ],
      version: 'v11',
    },
  };

  await zenodoRequest('PUT', '/api/records/' + draftId + '/draft', JSON.stringify(metadata), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const published = await zenodoRequest('POST', '/api/records/' + draftId + '/draft/actions/publish', '', 'application/json');
  console.log('  PUBLISHED!');
  console.log('  DOI: ' + (published.doi || 'check zenodo'));
  console.log('  URL: ' + (published.links?.self_html || published.links?.html || 'check zenodo'));
  console.log('\n=== Zenodo v11 upload complete ===');
}

main().catch(e => {
  console.error('Fatal: ' + e.message);
  process.exit(1);
});
