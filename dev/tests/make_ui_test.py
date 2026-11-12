import os

os.makedirs('dev/tmp', exist_ok=True)
with open('dev/tmp/test_default.html', 'rb') as f:
    html = f.read()
with open('dev/tests/ui_test.js', 'rb') as f:
    js = f.read()
probe = b'<script>\n' + js + b'\n</script>'
html = html.replace(b'</body>', probe + b'</body>', 1)
with open('dev/tmp/ui_test2.html', 'wb') as f:
    f.write(html)
print('wrote dev/tmp/ui_test2.html', len(html), 'bytes')
