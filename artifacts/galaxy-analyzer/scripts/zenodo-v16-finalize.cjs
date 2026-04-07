#!/usr/bin/env node
const fs = require('fs');
const https = require('https');
const pathMod = require('path');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const DRAFT_ID = '19446498';
const ARCHIVE = pathMod.join(__dirname, '..', 'zenodo-archives', 'galaxy-rotation-curve-v16.tar.gz');
const BASE = 'zenodo.org';
const UA = 'GalaxyRotationCurveAnalyzer/1.0 (mailto:antminers99@gmail.com)';

function req(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE, path: urlPath, method,
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Accept': 'application/json',
        'User-Agent': UA,
      },
    };
    if (body && contentType) {
      opts.headers['Content-Type'] = contentType;
      if (typeof body === 'string') opts.headers['Content-Length'] = Buffer.byteLength(body);
    }
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('  HTTP ' + res.statusCode + ': ' + d.substring(0, 500));
          reject(new Error('HTTP ' + res.statusCode));
        } else {
          try { resolve(JSON.parse(d)); } catch { resolve(d); }
        }
      });
    });
    r.on('error', reject);
    if (body && typeof body !== 'string') body.pipe(r);
    else { if (body) r.write(body); r.end(); }
  });
}

function uploadToBucket(bucketUrl, localPath, remoteName) {
  return new Promise((resolve, reject) => {
    const stat = fs.statSync(localPath);
    const url = new URL(bucketUrl + '/' + remoteName);
    const opts = {
      hostname: url.hostname,
      path: url.pathname,
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
        'User-Agent': UA,
      },
    };
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) {
          console.error('  Upload HTTP ' + res.statusCode + ': ' + d.substring(0, 300));
          reject(new Error('Upload fail'));
        } else {
          console.log('  Uploaded: ' + remoteName + ' (' + (stat.size / 1024 / 1024).toFixed(1) + ' MB)');
          resolve();
        }
      });
    });
    r.on('error', reject);
    fs.createReadStream(localPath).pipe(r);
  });
}

async function main() {
  console.log('=== Zenodo v16 Finalize ===\n');

  if (!fs.existsSync(ARCHIVE)) {
    console.error('Archive not found: ' + ARCHIVE);
    process.exit(1);
  }
  console.log('Archive: ' + ARCHIVE);
  console.log('Size: ' + (fs.statSync(ARCHIVE).size / 1024 / 1024).toFixed(1) + ' MB\n');

  console.log('Step 1: Fetch draft...');
  const draft = await req('GET', '/api/deposit/depositions/' + DRAFT_ID, null, null);
  console.log('  State: ' + draft.state);
  console.log('  Title: ' + draft.title);
  const bucketUrl = draft.links?.bucket;
  console.log('  Bucket: ' + bucketUrl);

  console.log('\nStep 2: Delete existing files...');
  const files = await req('GET', '/api/deposit/depositions/' + DRAFT_ID + '/files', null, null);
  console.log('  Found ' + files.length + ' existing files');
  for (const f of files) {
    console.log('    Deleting: ' + f.filename);
    await req('DELETE', '/api/deposit/depositions/' + DRAFT_ID + '/files/' + f.id, null, null);
  }

  console.log('\nStep 3: Upload archive via bucket...');
  await uploadToBucket(bucketUrl, ARCHIVE, 'galaxy-rotation-curve-v16.tar.gz');

  console.log('\nStep 4: Update metadata...');
  const meta = {
    metadata: {
      title: 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v16 — Programs 1-12, Dark Matter Discrimination)',
      upload_type: 'dataset',
      publication_date: draft.metadata?.publication_date || '2026-04-07',
      description:
        '<p><b>v16 — Dark matter model discrimination (Program 12: DM-1 through DM-4V).</b></p>' +
        '<p>Single archive (3.2 MB, 648 files): analysis scripts, result JSONs, source code, data files, documentation.</p>' +
        '<p><b>Core Result (unchanged):</b> Hidden state H drives bilateral Vf_Resid-a0_Resid coupling at r = 0.77 (LOO, p &lt; 0.001, N = 55), replicated on N = 59 independent galaxies. 1D information ceiling: 88.5% inaccessible. Red team: 11/11 pass.</p>' +
        '<p><b>Program 12 — Dark Matter Discrimination Series:</b></p>' +
        '<ul>' +
        '<li><b>DM-1 Kill Constraints:</b> 12 mandatory constraints tested against 7 DM models. CDM smooth DEAD (6 fails), MOND DEAD (10 fails). CDM+halo shape: 12/12 pass.</li>' +
        '<li><b>DM-2V Verification:</b> CDM+shape vs SIDM — SIDM signal was coverage artifact. r(DQ, outer C) = 0.746 vs r(DQ, inner C) = 0.431. LOO inner-led: 0/7. CDM+shape CONFIRMED.</li>' +
        '<li><b>DM-3 Fuzzy DM:</b> Wave fingerprint test — 0/4 for Fuzzy DM. No wave signature. Fuzzy DM DEAD.</li>' +
        '<li><b>DM-4 + DM-4V:</b> shapeAmplitude captures H at r = 0.80, p = 0.028, LOO 7/7 positive, bootstrap 97.8% positive, bar-exclusion stable. CONFIRMED.</li>' +
        '</ul>' +
        '<p><b>Model Status:</b> DEAD: CDM smooth, MOND, Fuzzy DM. No signal: SIDM, WDM. LEADING: CDM + non-axisymmetric halo shape (shapeAmplitude r=0.80).</p>' +
        '<p>Data: SPARC (Lelli, McGaugh &amp; Schombert 2016); THINGS (Walter et al. 2008).</p>',
      creators: draft.metadata?.creators || [{ name: 'Fnd89', affiliation: 'Independent Researcher' }],
      access_right: draft.metadata?.access_right || 'open',
      license: draft.metadata?.license || 'cc-by-4.0',
      version: 'v16',
      keywords: draft.metadata?.keywords || ['galaxy rotation curves', 'dark matter', 'SPARC', 'THINGS', 'radial acceleration relation'],
    },
  };
  await req('PUT', '/api/deposit/depositions/' + DRAFT_ID, JSON.stringify(meta), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const pub = await req('POST', '/api/deposit/depositions/' + DRAFT_ID + '/actions/publish', '', 'application/json');
  console.log('\n=== PUBLISHED ===');
  console.log('DOI: ' + (pub.doi || 'check zenodo'));
  console.log('URL: ' + (pub.links?.html || pub.links?.record_html || 'check zenodo'));
}

main().catch(e => { console.error('\nFAILED: ' + e.message); process.exit(1); });
