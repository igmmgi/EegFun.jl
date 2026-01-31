import { defineConfig } from 'vitepress'

export default defineConfig({
    ignoreDeadLinks: true,  // Allow building with TODO/placeholder links
    title: 'EegFun.jl',
    description: 'High-performance EEG data analysis in Julia',

    themeConfig: {
        nav: [
            { text: 'Home', link: '/' },
            { text: 'Tutorials', link: '/tutorials/getting-started' },
            { text: 'Reference', link: '/reference/index' },
            { text: 'GitHub', link: 'https://github.com/igmmgi/EegFun.jl' }
        ],

        sidebar: [
            {
                text: 'Tutorials',
                collapsed: false,
                items: [
                    { text: 'Getting Started', link: '/tutorials/getting-started' },
                    { text: 'Basic Preprocessing', link: '/tutorials/basic-preprocessing' },
                    { text: 'ERP Analysis', link: '/tutorials/erp-analysis' },
                    { text: 'ICA Workflow', link: '/tutorials/ica-workflow' },
                    { text: 'Statistical Analysis', link: '/tutorials/statistical-analysis' }
                ]
            },
            {
                text: 'How-To Guides',
                collapsed: false,
                items: [
                    { text: 'Filter Data', link: '/how-to/filter-data' },
                    { text: 'Create Epochs', link: '/how-to/create-epochs' },
                    { text: 'Topographic Plots', link: '/how-to/topographic-plots' }
                ]
            },
            {
                text: 'Explanations',
                collapsed: false,
                items: [
                    { text: 'Data Structures', link: '/explanations/data-structures' },
                    { text: 'ICA', link: '/explanations/ica' },
                    { text: 'Statistics', link: '/explanations/statistics' },
                    { text: 'Filtering', link: '/explanations/filtering' },
                    { text: 'ERP Measurements', link: '/explanations/erp-measurements' },
                    { text: 'Visualization', link: '/explanations/visualization' }
                ]
            },
            {
                text: 'API Reference',
                collapsed: false,
                items: [
                    { text: 'Overview', link: '/reference/index' },
                    { text: 'Data Loading', link: '/reference/data-loading' },
                    { text: 'Preprocessing', link: '/reference/preprocessing' },
                    { text: 'Types', link: '/reference/types' },
                    { text: 'Epochs', link: '/reference/epochs' },
                    { text: 'ERP', link: '/reference/erp' },
                    { text: 'ICA', link: '/reference/ica' },
                    { text: 'Statistics', link: '/reference/statistics' },
                    { text: 'Plotting', link: '/reference/plotting' },
                    { text: 'Layouts', link: '/reference/layouts' },
                    { text: 'Patterns', link: '/reference/patterns' }
                ]
            }
        ],

        socialLinks: [
            { icon: 'github', link: 'https://github.com/igmmgi/EegFun.jl' }
        ],

        search: {
            provider: 'local'
        },

        footer: {
            message: 'Released under the MIT License.',
            copyright: 'Copyright © 2024-present'
        }
    }
})
