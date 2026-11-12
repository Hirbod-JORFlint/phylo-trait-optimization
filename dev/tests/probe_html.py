import sys

src, dst = sys.argv[1], sys.argv[2]
with open(src, 'rb') as f:
    html = f.read()

if b'window.__errs' not in html:
    head_probe = (b"<script>window.__errs=[];window.addEventListener('error',function(e){"
                  b"window.__errs.push(String(e.message||e.error));});</script>")
    html = html.replace(b'</head>', head_probe + b'</head>', 1)

probe = b'''<script>
setTimeout(function(){
  var el=document.querySelector('.js-plotly-plot');
  var w=document.querySelector('.html-widget');
  var wc=document.getElementById('htmlwidget_container');
  function info(n){ if(!n) return null; var cs=getComputedStyle(n); return {tag:n.tagName,id:n.id,cls:n.className,style:n.getAttribute('style'),clientH:n.clientHeight,offsetH:n.offsetHeight,cssH:cs.height}; }
  var out={errs:window.__errs, plotly:typeof Plotly, hw:typeof HTMLWidgets};
  out.el=info(w); out.container=info(wc);
  out.bodyH=document.body.clientHeight; out.viewportH=window.innerHeight;
  if(el){var r=el.getBoundingClientRect();
    out.rect={w:Math.round(r.width),h:Math.round(r.height)};
    out.client={w:el.clientWidth,h:el.clientHeight};
    out.full=el._fullLayout?{w:el._fullLayout.width,h:el._fullLayout.height}:null;}
  out.uis=document.querySelectorAll('.itv-ui').length;
  out.widgets=document.querySelectorAll('.html-widget').length;
  out.plots=document.querySelectorAll('.js-plotly-plot').length;
  document.body.setAttribute('data-debug', JSON.stringify(out));
},6000);
</script>'''
html = html.replace(b'</body>', probe + b'</body>', 1)
with open(dst, 'wb') as f:
    f.write(html)
print('wrote', dst, len(html), 'bytes')
