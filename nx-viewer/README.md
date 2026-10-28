# NX CAD Viewer

Static GitHub Pages 3D viewer for STL and GLB files. Works on mobile and desktop.

**Native Siemens NX PRT files are not supported yet.** A licensed or otherwise compatible NX geometry reader must be integrated; simply renaming a PRT file will not work.

The site processes supported model files in the browser without uploading them to a server.

Deploy using the workflow in .github/workflows/nx-viewer-pages.yml and configure GitHub Pages to use GitHub Actions. The workflow publishes the nx-viewer directory at the repository's Pages root.

The viewer uses Three.js from a CDN and requires an internet connection.