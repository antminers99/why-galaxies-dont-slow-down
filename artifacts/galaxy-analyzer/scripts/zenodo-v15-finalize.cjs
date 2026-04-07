#!/usr/bin/env node
const https = require('https');
const TOKEN = process.env.ZENODO_TOKEN;
if (!TOKEN) { console.error('ZENODO_TOKEN not set'); process.exit(1); }
const DRAFT_ID = '19446498';

function req(method, urlPath, body) {
  return new Promise((resolve, reject) => {
    const opts = { hostname: 'zenodo.org', path: urlPath, method, headers: { 'Authorization': 'Bearer ' + TOKEN, 'Accept': 'application/json' } };
    if (body) { opts.headers['Content-Type'] = 'application/json'; opts.headers['Content-Length'] = Buffer.byteLength(body); }
    const r = https.request(opts, res => { let d = ''; res.on('data', c => d += c); res.on('end', () => {
      if (res.statusCode >= 400) { console.error('HTTP ' + res.statusCode + ':', d.substring(0, 500)); reject(new Error('HTTP ' + res.statusCode)); }
      else { try { resolve(JSON.parse(d)); } catch { resolve(d); } }
    }); });
    r.on('error', reject); if (body) r.write(body); r.end();
  });
}

async function main() {
  console.log('Step 1: Fetch draft ' + DRAFT_ID + '...');
  const draft = await req('GET', '/api/records/' + DRAFT_ID + '/draft');
  console.log('  Current title:', draft.metadata?.title?.substring(0, 80));
  console.log('  Files:', draft.files?.count || '?');

  console.log('Step 2: Update metadata...');
  const meta = draft.metadata;
  meta.title = 'Per-Galaxy a0 Variation: Angular Velocity-Field Complexity as Carrier of Hidden State (v15 - Programs 1-11, post-submission verification)';
  meta.version = 'v15';
  meta.description = '<p><b>v15 - Post-submission interpretive refinement.</b> This version preserves the core hidden-state result but updates the interpretation of the 2D carrier from a uniquely isolated m=2 mode to a broader angular/non-axisymmetric velocity-field complexity after Program 11 verification. Test A: radial m=2 gradient is universal. Test B: m=3 dominates non-rotational power (58% vs 11% for m=2). H tracks total gravitational complexity. Core result UNCHANGED. See CHANGELOG-v15.md.</p>';
  await req('PUT', '/api/records/' + DRAFT_ID + '/draft', JSON.stringify({ metadata: meta }));
  console.log('  Metadata updated.');

  console.log('Step 3: Publish...');
  const pub = await req('POST', '/api/records/' + DRAFT_ID + '/draft/actions/publish', '');
  console.log('=== PUBLISHED ===');
  console.log('DOI:', pub.doi);
  console.log('URL:', pub.links?.self_html || pub.links?.html);
}
main().catch(e => { console.error('FAILED:', e.message); process.exit(1); });
