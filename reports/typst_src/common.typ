#set page(paper: "a4", margin: (left: 25mm, right: 25mm, top: 24mm, bottom: 24mm), numbering: "第 1 页")
#set text(lang: "zh", font: ("Source Han Serif SC", "思源宋體"), size: 12pt, weight: 400)
#set par(justify: true, first-line-indent: (amount: 2em, all: true), leading: 0.86em, spacing: 0.58em)
#set heading(numbering: "1.")

#let mono = ("Consolas", "Source Han Sans SC", "思源黑體")
#let pale = rgb("f8f5ee")
#let line = rgb("cfc4b2")

#let cover(title, subtitle) = [
  #align(center)[
    #text(font: ("Source Han Sans SC", "思源黑體"), size: 19pt, weight: 700)[《控制系统仿真与 CAD》课程设计小报告]

    #v(0.9cm)
    #text(font: ("Source Han Sans SC", "思源黑體"), size: 17pt, weight: 700)[#title]

    #v(0.5cm)
    #text(size: 12pt)[#subtitle]
  ]
  #v(0.8cm)
]

#let codefile(path) = block(width: 100%, fill: pale, stroke: 0.5pt + line, inset: 7pt, radius: 3pt, above: 0.7em, below: 0.9em, breakable: true)[
  #set par(first-line-indent: 0pt)
  #text(font: mono, size: 9pt, weight: 500)[#raw(read(path), lang: "matlab")]
]

#let fig(path, caption-text, width: 100%) = figure(
  image(path, width: width),
  caption: caption-text,
)
