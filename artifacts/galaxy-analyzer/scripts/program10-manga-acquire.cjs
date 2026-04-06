const fs = require('fs');
const path = require('path');
const https = require('https');
const http = require('http');

const outDir = path.join(__dirname, '..', 'data', 'multi-survey-2d', 'manga');
fs.mkdirSync(outDir, { recursive: true });

const targets = [
  { name: 'UGC06786', ra: 176.8258, dec: 60.1131 },
  { name: 'UGC09037', ra: 212.8263, dec: 13.6511 },
  { name: 'NGC4559', ra: 188.9908, dec: 27.9597 },
  { name: 'UGC02953', ra: 60.4042, dec: 35.3478 },
  { name: 'UGC03546', ra: 101.8738, dec: 65.5917 },
  { name: 'UGC03580', ra: 103.4896, dec: 70.3672 },
];

function httpGet(url) {
  return new Promise((resolve, reject) => {
    const mod = url.startsWith('https') ? https : http;
    mod.get(url, { timeout: 15000 }, (res) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        return httpGet(res.headers.location).then(resolve).catch(reject);
      }
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => resolve({ status: res.statusCode, data }));
    }).on('error', reject).on('timeout', function() { this.destroy(); reject(new Error('timeout')); });
  });
}

function httpGetBinary(url) {
  return new Promise((resolve, reject) => {
    const mod = url.startsWith('https') ? https : http;
    mod.get(url, { timeout: 30000 }, (res) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        return httpGetBinary(res.headers.location).then(resolve).catch(reject);
      }
      const chunks = [];
      res.on('data', chunk => chunks.push(chunk));
      res.on('end', () => resolve({ status: res.statusCode, data: Buffer.concat(chunks) }));
    }).on('error', reject).on('timeout', function() { this.destroy(); reject(new Error('timeout')); });
  });
}

async function findMaNGAPlateIFU(ra, dec, name) {
  const radius = 0.02;
  const url = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=' +
    encodeURIComponent(
      "SELECT TOP 1 p.mangaid, p.plate, p.ifudsgn, p.objra, p.objdec, p.nsa_elpetro_mass, p.nsa_z " +
      "FROM mangaDrpAll AS p " +
      "WHERE p.objra BETWEEN " + (ra - radius) + " AND " + (ra + radius) +
      " AND p.objdec BETWEEN " + (dec - radius) + " AND " + (dec + radius)
    ) + '&format=json';

  console.log('  Querying SDSS for ' + name + ' (RA=' + ra.toFixed(4) + ', Dec=' + dec.toFixed(4) + ')...');
  try {
    const resp = await httpGet(url);
    if (resp.status !== 200) {
      console.log('    HTTP ' + resp.status);
      return null;
    }
    const parsed = JSON.parse(resp.data);
    const rows = parsed[0]?.Rows || parsed;
    if (!rows || rows.length === 0) {
      console.log('    No MaNGA observation found within ' + (radius * 3600).toFixed(0) + '" radius');
      return null;
    }
    const row = rows[0];
    const mangaid = row.mangaid || row.MANGAID;
    const plate = row.plate || row.PLATE;
    const ifu = row.ifudsgn || row.IFUDSGN;
    console.log('    Found: mangaid=' + mangaid + '  plate=' + plate + '  ifu=' + ifu);
    return { mangaid, plate, ifu, ra: row.objra || row.OBJRA, dec: row.objdec || row.OBJDEC };
  } catch(e) {
    console.log('    Error: ' + e.message);
    return null;
  }
}

async function downloadMaNGAMap(plate, ifu, name) {
  const dapType = 'HYB10-MILESHC-MASTARSSP';
  const url = 'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/' +
    plate + '/' + ifu + '/manga-' + plate + '-' + ifu + '-MAPS-' + dapType + '.fits.gz';

  console.log('    Downloading DAP map: manga-' + plate + '-' + ifu + '...');
  try {
    const resp = await httpGetBinary(url);
    if (resp.status !== 200) {
      console.log('    HTTP ' + resp.status + ' — file not found or access denied');
      return null;
    }
    const outFile = path.join(outDir, name + '_MANGA_VFIELD.fits.gz');
    fs.writeFileSync(outFile, resp.data);
    console.log('    Saved: ' + outFile + ' (' + (resp.data.length / 1024 / 1024).toFixed(1) + ' MB)');
    return outFile;
  } catch(e) {
    console.log('    Download error: ' + e.message);
    return null;
  }
}


async function main() {
  console.log('='.repeat(72));
  console.log('PROGRAM 10 — MaNGA DATA ACQUISITION');
  console.log('='.repeat(72));

  const found = [];
  const notFound = [];

  for (const t of targets) {
    console.log('\n--- ' + t.name + ' ---');
    const info = await findMaNGAPlateIFU(t.ra, t.dec, t.name);
    if (info) {
      found.push({ ...t, ...info });
    } else {
      notFound.push(t);
    }
  }

  console.log('\n\n' + '#'.repeat(72));
  console.log('MaNGA CROSS-MATCH RESULTS');
  console.log('#'.repeat(72));
  console.log('\n  Found in MaNGA: ' + found.length + '/' + targets.length);
  for (const f of found) {
    console.log('    ' + f.name.padEnd(15) + 'mangaid=' + f.mangaid + '  plate=' + f.plate + '  ifu=' + f.ifu);
  }
  console.log('\n  Not found: ' + notFound.length);
  for (const nf of notFound) {
    console.log('    ' + nf.name);
  }

  if (found.length > 0) {
    console.log('\n\nAttempting to download DAP velocity maps...\n');
    for (const f of found) {
      await downloadMaNGAMap(f.plate, f.ifu, f.name);
    }
  }

  const catalogFile = path.join(outDir, 'manga-crossmatch.json');
  fs.writeFileSync(catalogFile, JSON.stringify({ found, notFound, timestamp: new Date().toISOString() }, null, 2));
  console.log('\n  Catalog saved to: ' + catalogFile);

  console.log('\n' + '='.repeat(72));
  console.log('ACQUISITION COMPLETE');
  console.log('='.repeat(72));
}

main().catch(e => console.error('Fatal:', e));
