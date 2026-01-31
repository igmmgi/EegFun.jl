import { defineConfig } from 'vitepress'
import path from 'path'
import { fileURLToPath } from 'url'

const __dirname = path.dirname(fileURLToPath(import.meta.url))

export default defineConfig({
    base: '/dev/1/',  // Critical: tells VitePress the base path for all links
    ignoreDeadLinks: true,

    vite: {
        resolve: {
            alias: {
                '@': path.resolve(__dirname, 'theme')
            }
        }
    },
})
