@echo off
rem Double-click to open the mechanism viewer in your browser.
rem Keep this window open while you use the viewer; close it when done.
cd /d "%~dp0"
start "" "http://localhost:8765/web/index.html"
".venv\Scripts\python.exe" -m http.server 8765
