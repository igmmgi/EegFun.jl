<script setup>
import { useData } from 'vitepress'
import { watchEffect } from 'vue'

const { page } = useData()
</script>

<template>
  <div v-if="page.headers" class="VPDocFooter">
      <div v-if="page.lastUpdated" class="last-updated">
        <p>Last updated: {{ new Date(page.lastUpdated).toLocaleDateString() }}</p>
      </div>
  </div>
</template>

<style scoped>
.VPDocFooter {
  margin-top: 64px;
}
</style>
