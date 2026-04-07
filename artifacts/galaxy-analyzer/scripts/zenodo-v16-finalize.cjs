#!/usr/bin/env node
const fs = require('fs');
const https = require('https');

const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }

const DRAFT_ID = '19446498';
const ARCHIVE = require('path').join(__dirname, '..', 'zenodo-archives', 'galaxy-rotation-curve-v16.tar.gz');
const BASE = 'zenodo.org';

function req(method, urlPath, body, contentType) {
  return new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE, path: urlPath, method,
      headers: { 'Authorization': 'Bearer ' + TOKEN, 'Accept': 'application/json' },
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

async function deleteFile(key) {
  await req('DELETE', '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(key), null, null);
}

async function uploadFile(localPath, remoteName) {
  const stat = fs.statSync(localPath);
  await req('POST', '/api/records/' + DRAFT_ID + '/draft/files',
    JSON.stringify([{ key: remoteName }]), 'application/json');
  await new Promise((resolve, reject) => {
    const opts = {
      hostname: BASE,
      path: '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(remoteName) + '/content',
      method: 'PUT',
      headers: {
        'Authorization': 'Bearer ' + TOKEN,
        'Content-Type': 'application/octet-stream',
        'Content-Length': stat.size,
      },
    };
    const r = https.request(opts, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => {
        if (res.statusCode >= 400) { console.error('Upload err:', res.statusCode, d.substring(0, 200)); reject(new Error('fail')); }
        else resolve();
      });
    });
    r.on('error', reject);
    fs.createReadStream(localPath).pipe(r);
  });
  await req('POST', '/api/records/' + DRAFT_ID + '/draft/files/' + encodeURIComponent(remoteName) + '/commit',
    '', 'application/json');
  console.log('  Uploaded: ' + remoteName + ' (' + (stat.size / 1024 / 1024).toFixed(1) + ' MB)');
}

async function main() {
  console.log('=== Zenodo v16 Finalize: Single Archive ===\n');
  console.log('Archive: ' + ARCHIVE);
  console.log('Draft ID: ' + DRAFT_ID + '\n');

  if (!fs.existsSync(ARCHIVE)) {
    console.error('Archive not found: ' + ARCHIVE);
    process.exit(1);
  }

  console.log('Step 1: Fetch draft...');
  const draft = await req('GET', '/api/records/' + DRAFT_ID + '/draft', null, null);
  console.log('  Title: ' + (draft.metadata?.title || 'unknown'));

  console.log('\nStep 2: Delete all existing files...');
  const files = await req('GET', '/api/records/' + DRAFT_ID + '/draft/files', null, null);
  const entries = files.entries || [];
  console.log('  Found ' + entries.length + ' files to delete');
  for (const f of entries) {
    console.log('    Deleting: ' + f.key);
    await deleteFile(f.key);
  }
  console.log('  Done.');

  console.log('\nStep 3: Upload archive...');
  await uploadFile(ARCHIVE, 'galaxy-rotation-curve-v16.tar.gz');

  console.log('\nStep 4: Update metadata...');
  const meta = draft.metadata;
  meta.title = 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v16 — Programs 1-12, Dark Matter Discrimination)';
  meta.version = 'v16';
  meta.description =
    '<p><b>v16 — Dark matter model discrimination (Program 12: DM-1 through DM-4V).</b></p>' +
    '<p>Single archive (3.2 MB, 646 files): analysis scripts, result JSONs, ' +
    'source code, data files, documentation.</p>' +
    '<p><b>Core Result (unchanged):</b> Hidden state H drives bilateral Vf_Resid-a0_Resid coupling ' +
    'at r = 0.77 (LOO, p &lt; 0.001, N = 55), replicated on N = 59 independent galaxies. ' +
    '1D information ceiling: 88.5% inaccessible. Red team: 11/11 pass.</p>' +
    '<p><b>Program 12 — Dark Matter Discrimination Series:</b></p>' +
    '<ul>' +
    '<li><b>DM-1 Kill Constraints:</b> 12 mandatory constraints tested against 7 DM models. ' +
    'CDM smooth DEAD (6 fails), MOND DEAD (10 fails). CDM+halo shape: 12/12 pass.</li>' +
    '<li><b>DM-2V Verification:</b> CDM+shape vs SIDM decisive test. ' +
    'SIDM signal was coverage artifact — r(DQ, outer C) = 0.746 vs r(DQ, inner C) = 0.431, ' +
    'LOO inner-led: 0/7, bootstrap inner-led: 19.7%. CDM+shape CONFIRMED.</li>' +
    '<li><b>DM-3 Fuzzy DM:</b> Wave fingerprint test — 0/4 for Fuzzy DM. ' +
    'No wave signature (entropy 0.72, flip rate 0.18, odd/even 4.15). Fuzzy DM DEAD.</li>' +
    '<li><b>DM-4 + DM-4V:</b> Quantitative halo shape — shapeAmplitude captures H at r = 0.80, ' +
    'p = 0.028, LOO 7/7 positive, bootstrap 97.8% positive, bar-exclusion stable. CONFIRMED.</li>' +
    '</ul>' +
    '<p><b>Model Status:</b> DEAD: CDM smooth, MOND, Fuzzy DM. ' +
    'No signal: SIDM, WDM. LEADING: CDM + non-axisymmetric halo shape (shapeAmplitude r=0.80).</p>' +
    '<p>Data: SPARC (Lelli, McGaugh &amp; Schombert 2016); THINGS (Walter et al. 2008).</p>';

  await req('PUT', '/api/records/' + DRAFT_ID + '/draft', JSON.stringify({ metadata: meta }), 'application/json');
  console.log('  Metadata updated.');

  console.log('\nStep 5: Publish...');
  const pub = await req('POST', '/api/records/' + DRAFT_ID + '/draft/actions/publish', '', 'application/json');
  console.log('\n=== PUBLISHED ===');
  console.log('DOI: ' + (pub.doi || 'check zenodo'));
  console.log('URL: ' + (pub.links?.self_html || pub.links?.html || 'check zenodo'));
}

main().catch(e => { console.error('\nFAILED: ' + e.message); process.exit(1); });
