<script setup>
defineProps({
  text: {
    type: String,
    required: true
  }
})
</script>

<template>
  <span class="badge">{{ text }}</span>
</template>

<style scoped>
.badge {
  display: inline-block;
  padding: 0 0.5em;
  border-radius: 4px;
  background-color: var(--vp-c-brand);
  color: white;
  font-size: 0.8em;
  font-weight: bold;
}
</style>
