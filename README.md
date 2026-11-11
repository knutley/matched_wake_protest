# Partisan Wakes

Companion site for "Does Partisan Bias Become Institutionalised?" (Katelyn Nutley, working paper).

- `index.html`: landing page
- `wake/index.html`: interactive matched wake analysis example (one UK pair, illustrative events)
- `talk/talk.qmd`: conference slides (Quarto reveal.js)
- `talk/figures/`: put `fig_1_main.png` and `fig_2_baselines.png` from the paper here

## Publishing

1. Create a new GitHub repository (e.g. `partisan-wakes`) and upload these files, including the hidden `.github` folder.
2. In the repository: **Settings → Pages → Build and deployment → Source: GitHub Actions**.
3. Every push to `main` re-renders the slides and republishes the site at
   `https://<username>.github.io/<repository>/`.

## Editing the slides locally

```
quarto preview talk/talk.qmd
```

PDF backup for conferences: open the slides in Chrome, press `E`, then print to PDF.

## Data

ACLED event data are not included and must not be committed. The worked example uses illustrative events.
