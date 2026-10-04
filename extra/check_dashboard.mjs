/*
 * Dashboard smoke test. Loads a built dashboard in a DOM, clicks every tab and control, and
 * audits the SVG geometry — a label that runs past its chart box is a real defect the palette
 * validator cannot catch.
 *
 *   npm install jsdom && node extra/check_dashboard.mjs results_nf/dashboard/index.html
 */
import { JSDOM, VirtualConsole } from 'jsdom';
import fs from 'fs';

const html = fs.readFileSync(process.argv[2], 'utf8');
const errors = [];
const vc = new VirtualConsole();
vc.on('jsdomError', e => errors.push('jsdomError: ' + (e.stack || e.message)));
vc.on('error', (...a) => errors.push('console.error: ' + a.join(' ')));
vc.on('warn', (...a) => errors.push('console.warn: ' + a.join(' ')));

const dom = new JSDOM(html, { runScripts: 'dangerously', pretendToBeVisual: true, virtualConsole: vc });
const { window } = dom;
const doc = window.document;

const tabs = [...doc.querySelectorAll('nav.tabs button')];
console.log('tabs:', tabs.map(t => t.textContent).join(' | '));

function summarize(id) {
  const p = doc.getElementById('panel-' + id);
  const cards = p.querySelectorAll('.card').length;
  const svgs = p.querySelectorAll('svg').length;
  const marks = p.querySelectorAll('svg path, svg circle').length;
  const tables = p.querySelectorAll('table').length;
  const rows = p.querySelectorAll('tbody tr').length;
  const legends = p.querySelectorAll('.legend').length;
  const heads = [...p.querySelectorAll('h2, .figcap h3')].map(x => x.textContent);
  const empty = [...p.querySelectorAll('.empty')].map(x => x.textContent);
  console.log(`\n[${id}] cards=${cards} svg=${svgs} marks=${marks} tables=${tables} tbodyRows=${rows} legends=${legends}`);
  heads.forEach(t => console.log('    · ' + t));
  empty.forEach(t => console.log('    EMPTY: ' + t));
  // geometry audit: no text may run past either edge of its 760-unit chart box
  const bad = [];
  p.querySelectorAll('svg').forEach(svg => {
    svg.querySelectorAll('text').forEach(t => {
      if (t.getAttribute('transform')) return;   // rotated axis titles are measured on the other axis
      const w = (t.textContent || '').length * 5.6;
      const x = parseFloat(t.getAttribute('x') || '0');
      const anchor = t.getAttribute('text-anchor');
      const left = anchor === 'end' ? x - w : anchor === 'middle' ? x - w / 2 : x;
      const right = left + w;
      if (left < -1 || right > 761) bad.push(`"${(t.textContent||'').slice(0,32)}" [${left.toFixed(0)},${right.toFixed(0)}]`);
    });
  });
  if (bad.length) { console.log('    TEXT OVERFLOW:', bad.length); bad.slice(0,5).forEach(b => console.log('      ' + b)); }
}

for (const t of tabs) { t.click(); summarize(t.id.replace('tab-', '')); }

// interactions
console.log('\n--- interactions ---');
doc.getElementById('tab-expression').click();
const selects = [...doc.getElementById('panel-expression').querySelectorAll('select')];
console.log('expression selects:', selects.length);
if (selects[0]) {
  const before = doc.getElementById('panel-expression').textContent.length;
  selects[0].value = selects[0].options[selects[0].options.length - 1].value;
  selects[0].dispatchEvent(new window.Event('change'));
  console.log('contrast switch ->', selects[0].value, 'panel text len', before, '→', doc.getElementById('panel-expression').textContent.length);
}
doc.getElementById('tab-genes').click();
const gp = doc.getElementById('panel-genes');
const search = gp.querySelector('input[type=search]');
if (search) {
  const before = gp.querySelectorAll('tbody tr').length;
  search.value = 'ribosom';
  search.dispatchEvent(new window.Event('input'));
  console.log('gene search "ribosom": rows', before, '→', gp.querySelectorAll('tbody tr').length);
  const th = [...gp.querySelectorAll('th')].find(x => x.textContent.startsWith('log2FC'));
  th.dispatchEvent(new window.MouseEvent('click', { bubbles: true }));
  const tbs = gp.querySelectorAll('tbody'); console.log('sorted by', th.textContent.trim(), '| first row:', tbs[tbs.length-1].querySelector('tr')?.textContent.slice(0, 70));
}
// table twin toggles
doc.getElementById('tab-qc').click();
const toggles = [...doc.getElementById('panel-qc').querySelectorAll('button.ghost')].filter(b => b.textContent === 'Table');
console.log('table toggles in QC:', toggles.length);
if (toggles[0]) { toggles[0].click(); console.log('after toggle, chart hidden =', !!doc.getElementById('panel-qc').querySelector('div[hidden] svg')); }
// theme
const tb = doc.getElementById('theme');
tb.click(); console.log('theme ->', doc.documentElement.getAttribute('data-theme'), '/', tb.textContent);
tb.click(); console.log('theme ->', doc.documentElement.getAttribute('data-theme'), '/', tb.textContent);
// methods rail
doc.getElementById('tab-overview').click();
const rail = [...doc.getElementById('panel-overview').querySelectorAll('.rail button')];
console.log('rail stages:', rail.length, '| first:', rail[0]?.textContent);
rail[3]?.click();
console.log('after rail click, open stage:', [...doc.querySelectorAll('details.stage[open]')].map(d => d.id).join(','));

console.log('\nERRORS:', errors.length);
errors.slice(0, 12).forEach(e => console.log('  ' + e.slice(0, 400)));
