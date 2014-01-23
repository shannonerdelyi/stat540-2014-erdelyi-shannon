Playing with Markdown
=====================
#### STAT 540: Seminar 3
#### Shannon Erdelyi
#### 22-01-2014

I always use *Sweave* to create documents with **R** because I like pdf files. I have been putting off learning markdown, even though Jenny always says:
> Markdown is not that hard to learn!

I can think of a few reasons why markdown is better than sweave:

* It's web friendly
* I can still write nice LaTeX equations
  * $\alpha=2$
  * $\frac{a}{b}$
* It's readable 

Let's try some r code.

```r
a <- 1
b <- 2
(f <- a/b)
```

```
## [1] 0.5
```


I can also plot stuff.
### This is my simple plot

```r
plot(a, b)
```

![plot of chunk simplePlot](figure/simplePlot.png) 


### My only complaint
I don't like all of the files created in my working directory when I knit2html. Can I get rid of these somehow?
