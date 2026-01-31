import { defineConfig } from 'vitepress'

export default defineConfig({
    ignoreDeadLinks: true,  // Allow building with TODO/placeholder links
    title: 'EegFun.jl',
    description: 'High-performance EEG data analysis in Julia',

    themeConfig: {
        logo: '/logo.png',  // Logo appears in navbar

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
                text: 'Demos',
                collapsed: true,
                items: [
                    {
                        text: 'Data Processing',
                        collapsed: true,
                        items: [
                            { text: 'Baseline Correction', link: '/demos/baseline' },
                            { text: 'Mirror Padding', link: '/demos/mirror' },
                            { text: 'Rereferencing', link: '/demos/rereference' },
                            { text: 'Resampling', link: '/demos/resample' },
                            { text: 'Data Loading', link: '/demos/data' }
                        ]
                    },
                    {
                        text: 'Channel Operations',
                        collapsed: true,
                        items: [
                            { text: 'Channel Metrics', link: '/demos/channel-metrics' },
                            { text: 'Channel Repair', link: '/demos/channel-repair' },
                            { text: 'Channel Summary', link: '/demos/channel-summary' }
                        ]
                    },
                    {
                        text: 'Artifact Detection',
                        collapsed: true,
                        items: [
                            { text: 'Artifacts', link: '/demos/artifacts' },
                            { text: 'Plot Artifacts', link: '/demos/plot-artifacts' },
                            { text: 'Joint Probability', link: '/demos/plot-joint-probability' }
                        ]
                    },
                    { text: 'ICA Analysis', link: '/demos/ica' },
                    {
                        text: 'Time-Frequency',
                        collapsed: true,
                        items: [
                            { text: 'Morlet Wavelets', link: '/demos/tf-morlet' },
                            { text: 'Multitaper', link: '/demos/tf-multitaper' },
                            { text: 'STFT', link: '/demos/tf-stft' }
                        ]
                    },
                    {
                        text: 'Statistical Analysis',
                        collapsed: true,
                        items: [
                            { text: 'ERP Measurements', link: '/demos/erp-measurements' },
                            { text: 'Statistics', link: '/demos/statistics' },
                            { text: 'Decoding', link: '/demos/decoding' },
                            { text: 'RSA', link: '/demos/rsa' }
                        ]
                    },
                    {
                        text: 'Visualization',
                        collapsed: true,
                        items: [
                            { text: 'Channel Spectrum', link: '/demos/plot-channel-spectrum' },
                            { text: 'Channel Summary', link: '/demos/plot-channel-summary' },
                            { text: 'Correlation Heatmap', link: '/demos/plot-correlation-heatmap' },
                            { text: 'Data Browser', link: '/demos/plot-databrowser' },
                            { text: 'Epochs', link: '/demos/plot-epochs' },
                            { text: 'ERP', link: '/demos/plot-erp' },
                            { text: 'ERP Image', link: '/demos/plot-erp-image' },
                            { text: 'Filter Response', link: '/demos/plot-filter' },
                            { text: 'Layout', link: '/demos/plot-layout' },
                            { text: 'Topography', link: '/demos/plot-topography' },
                            { text: 'Triggers', link: '/demos/plot-triggers' }
                        ]
                    }
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
