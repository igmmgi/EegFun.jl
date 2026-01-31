import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/igmmgi.github.io/EegFun.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'Home', link: '/index' },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Getting Started', link: '/tutorials/getting-started' },
{ text: 'Basic Preprocessing', link: '/tutorials/basic-preprocessing' },
{ text: 'ERP Analysis', link: '/tutorials/erp-analysis' },
{ text: 'ICA Workflow', link: '/tutorials/ica-workflow' }]
 },
{ text: 'How-To Guides', collapsed: false, items: [
{ text: 'Filter Data', link: '/how-to/filter-data' },
{ text: 'Create Epochs', link: '/how-to/create-epochs' },
{ text: 'Topographic Plots', link: '/how-to/topographic-plots' }]
 },
{ text: 'Explanations', collapsed: false, items: [
{ text: 'Data Structures', link: '/explanations/data-structures' },
{ text: 'Statistical Methods', link: '/explanations/statistics' }]
 },
{ text: 'Reference', collapsed: false, items: [
{ text: 'Overview', link: '/reference/index' },
{ text: 'Types', link: '/reference/types' }]
 }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/igmmgi.github.io/EegFun.jl/dev/',// TODO: replace this in makedocs!
  title: 'EegFun.jl',
  description: 'Documentation for EegFun.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('/igmmgi.github.io/EegFun.jl'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'Home', link: '/index' },
{ text: 'Tutorials', collapsed: false, items: [
{ text: 'Getting Started', link: '/tutorials/getting-started' },
{ text: 'Basic Preprocessing', link: '/tutorials/basic-preprocessing' },
{ text: 'ERP Analysis', link: '/tutorials/erp-analysis' },
{ text: 'ICA Workflow', link: '/tutorials/ica-workflow' }]
 },
{ text: 'How-To Guides', collapsed: false, items: [
{ text: 'Filter Data', link: '/how-to/filter-data' },
{ text: 'Create Epochs', link: '/how-to/create-epochs' },
{ text: 'Topographic Plots', link: '/how-to/topographic-plots' }]
 },
{ text: 'Explanations', collapsed: false, items: [
{ text: 'Data Structures', link: '/explanations/data-structures' },
{ text: 'Statistical Methods', link: '/explanations/statistics' }]
 },
{ text: 'Reference', collapsed: false, items: [
{ text: 'Overview', link: '/reference/index' },
{ text: 'Types', link: '/reference/types' }]
 }
]
,
    editLink: { pattern: "https://https://github.com/igmmgi/EegFun.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/igmmgi/EegFun.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
