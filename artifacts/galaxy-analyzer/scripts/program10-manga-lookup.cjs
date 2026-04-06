const https = require('https');
const fs = require('fs');
const path = require('path');

function httpGet(url) {
  return new Promise((resolve, reject) => {
    https.get(url, { timeout: 20000 }, (res) => {
      if (res.statusCode === 301 || res.statusCode === 302) {
        return httpGet(res.headers.location).then(resolve).catch(reject);
      }
      let data = '';
      res.on('data', chunk => data += chunk);
      res.on('end', () => resolve({ status: res.statusCode, data }));
    }).on('error', reject).on('timeout', function() { this.destroy(); reject(new Error('timeout')); });
  });
}

const targets = [
  { name: 'NGC2841', ra: 140.5108, dec: 50.9764 },
  { name: 'UGC06786', ra: 176.826, dec: 60.113 },
  { name: 'UGC09037', ra: 212.141, dec: 13.694 },
  { name: 'NGC4559', ra: 188.990, dec: 27.960 },
  { name: 'UGC02953', ra: 60.770, dec: 35.350 },
  { name: 'NGC5055', ra: 198.956, dec: 42.029 },
  { name: 'NGC3521', ra: 166.452, dec: -0.036 },
  { name: 'NGC7331', ra: 339.267, dec: 34.416 },
  { name: 'NGC2403', ra: 114.214, dec: 65.603 },
  { name: 'NGC3198', ra: 154.979, dec: 45.550 },
  { name: 'UGC03546', ra: 101.874, dec: 65.592 },
  { name: 'UGC03580', ra: 103.490, dec: 70.367 },
  { name: 'NGC4826', ra: 194.182, dec: 21.683 },
  { name: 'IC2574', ra: 157.089, dec: 68.413 },
  { name: 'NGC4088', ra: 182.576, dec: 50.540 },
];

async function querySDSS(ra, dec, name, radius) {
  const sql = "SELECT TOP 3 p.mangaid, p.plate, p.ifudsgn, p.objra, p.objdec, p.nsa_elpetro_mass, p.nsa_z " +
    "FROM mangaDrpAll AS p " +
    "WHERE p.objra BETWEEN " + (ra - radius) + " AND " + (ra + radius) +
    " AND p.objdec BETWEEN " + (dec - radius) + " AND " + (dec + radius);

  const url = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=' +
    encodeURIComponent(sql) + '&format=json';

  try {
    const resp = await httpGet(url);
    if (resp.status !== 200) return [];
    const parsed = JSON.parse(resp.data);
    if (Array.isArray(parsed) && parsed.length > 0 && parsed[0].Rows) return parsed[0].Rows;
    if (Array.isArray(parsed)) return parsed;
    return [];
  } catch(e) {
    return [];
  }
}

async function main() {
  console.log('='.repeat(72));
  console.log('PROGRAM 10 — MaNGA PLATE-IFU LOOKUP');
  console.log('Searching SDSS DR17 with wider radius (0.1 deg)');
  console.log('='.repeat(72));

  const found = [];
  for (const t of targets) {
    process.stdout.write('\n  ' + t.name.padEnd(12));
    let rows = await querySDSS(t.ra, t.dec, t.name, 0.1);
    if (rows.length === 0) {
      rows = await querySDSS(t.ra, t.dec, t.name, 0.2);
    }
    if (rows.length > 0) {
      const r = rows[0];
      const mangaid = r.mangaid || r.MANGAID || 'unknown';
      const plate = r.plate || r.PLATE || 0;
      const ifu = r.ifudsgn || r.IFUDSGN || 0;
      const z = r.nsa_z || r.NSA_Z || 0;
      console.log(' FOUND  mangaid=' + mangaid + '  plate=' + plate + '  ifu=' + ifu + '  z=' + (z || 0).toFixed ? z.toFixed(4) : z);
      found.push({ name: t.name, mangaid, plate, ifu, z });
    } else {
      console.log(' NOT IN MaNGA');
    }
    await new Promise(r => setTimeout(r, 500));
  }

  console.log('\n\n  SUMMARY: Found ' + found.length + '/' + targets.length + ' in MaNGA DR17');

  if (found.length > 0) {
    console.log('\n  MaNGA DAP URLs for download:');
    for (const f of found) {
      const url = 'https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/' +
        f.plate + '/' + f.ifu + '/manga-' + f.plate + '-' + f.ifu + '-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz';
      console.log('    ' + f.name + ': ' + url);
    }
  }

  fs.writeFileSync(
    path.join(__dirname, '..', 'data', 'multi-survey-2d', 'manga', 'manga-lookup.json'),
    JSON.stringify({ found, timestamp: new Date().toISOString() }, null, 2)
  );
}

main().catch(e => console.error('Fatal:', e));
