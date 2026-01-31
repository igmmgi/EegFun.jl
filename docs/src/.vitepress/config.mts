import { defineConfig } from 'vitepress'

export default defineConfig({
    base: '/dev/1/',  // Critical: tells VitePress the base path for all links
    ignoreDeadLinks: true,
})
