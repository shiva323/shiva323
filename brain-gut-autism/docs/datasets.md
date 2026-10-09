# دیتاست‌های محور روده–مغز (Gut–Brain Axis) در اوتیسم کودکان

> آخرین بررسی: اکتبر ۲۰۲۶. لینک‌ها و شرایط دسترسی را قبل از شروع دوباره چک کن؛
> مخزن‌ها و شرایط اشتراک‌گذاری داده زیاد تغییر می‌کنند.

## خلاصهٔ مهم (اول این را بخوان)

1. **هیچ دیتاست کاملاً عمومی و قابل‌دانلودی پیدا نکردم که برای *همان کودکان* اوتیسم هم
   EEG داشته باشد و هم میکروبیوم.** برای MRI + میکروبیوم هم دیتاست‌های خوبی منتشر
   شده‌اند، ولی داده‌شان باید **از نویسندگان درخواست** شود (دسترسی کنترل‌شده).
2. بنابراین سه مسیر داری (جزئیات در انتهای همین فایل):
   - **مسیر A:** درخواست دادهٔ «زوج» (paired) از گروه‌های Mi و همکاران (۲۰۲۶) و Aziz-Zadeh و همکاران (۲۰۲۵)؛
   - **مسیر B:** هم‌اکنون روی داده‌های باز کار کنی: شاخهٔ روده (میکروبیوم + آزمایش + مقیاس‌های بالینی)
     و شاخهٔ مغز (EEG/fMRI) را *جداگانه* بسازی و اعتبارسنجی کنی؛
   - **مسیر C:** داده‌ی خودت را جمع کنی (EEG استراحت + نمونهٔ مدفوع + پرسشنامه‌ها). با تخصص
     EEG تو، این عملی‌ترین راه برای یک مقالهٔ واقعاً «روده–مغز» است.
   - **مسیر D (بدون درخواست از نویسنده):** چهار طرح فقط با داده‌های کاملاً باز؛ بخش بعد را ببین.
3. **هشدار روش‌شناختی:** دو دیتاست که افرادشان متفاوت‌اند (مثلاً ABIDE برای مغز و یک کوهورت
   میکروبیوم دیگر) را **نمی‌توان در سطح فرد ادغام کرد**. ادغام (fusion) فقط وقتی معنی دارد
   که هر دو داده از همان کودک باشد.

---

## سطح ۰ — فقط داده‌های کاملاً باز (بدون نیاز به درخواست از نویسنده)

هیچ دیتاست بازی پیدا نکردم که برای **همان کودکان اوتیسم** هم مغز داشته باشد و هم روده. اما
چهار طرح زیر را می‌شود فقط با دانلود مستقیم انجام داد. در هر طرح، نوع ارتباطی که با اوتیسم
برقرار می‌شود جداگانه نوشته شده است.

### D1 — RESONANCE: میکروبیوم + MRI + آزمون‌های شناختی در ۳۸۱ کودک ⭐ پیشنهاد اول

| | |
|---|---|
| **مقاله** | Bonham et al., 2023, *Science Advances*، با عنوان *Gut-resident microorganisms and their genes are associated with cognition and neuroanatomy in children* |
| **افراد** | ۳۸۱ کودک با رشد معمول، از ۴۰ روزگی تا ۱۰ سالگی |
| **داده** | متاژنوم شاتگان مدفوع (پروفایل گونه‌ها، gene families، KO و Pfam)، **حجم نواحی مغز از MRI** (`brain_normalized.csv`)، نمرات شناختی و متادیتا (`complete_filtered_dataset.csv`) |
| **دسترسی** | دانلود آزاد از Dryad، حدود ۶.۲ گیگابایت: [doi:10.5061/dryad.fxpnvx0zq](https://datadryad.org/dataset/doi:10.5061/dryad.fxpnvx0zq) |

**ارتباط با اوتیسم: طرح «امضای میکروبی اوتیسم»**
1. با داده‌های باز اوتیسم (پایین‌تر آمده‌اند) یک امضای میکروبی ASD در برابر TD بساز. منظور مجموعه‌ای از گونه‌ها یا مسیرهای عملکردی است، مثل Bifidobacterium، Desulfovibrio، تولید SCFA، مسیر تریپتوفان و GABA.
2. به هر کودک RESONANCE یک «نمرهٔ شباهت به اوتیسم» بده.
3. بررسی کن آیا این نمره با حجم نواحی مرتبط با اوتیسم رابطه دارد، مثل اینسولا، سینگولیت، آمیگدال و مخچه. این کار را با کنترل سن، جنس و عمق توالی انجام بده. چون نمره و MRI از **همان کودک** هستند، این تحلیل در سطح فرد معتبر است.
4. محدودیت: این کودکان اوتیسم ندارند، پس نتیجه یک تحلیل **بُعدی (dimensional)** است و تشخیصی نیست.

پایپ‌لاین `gutbrain` بدون تغییر روی این داده کار می‌کند. پروفایل گونه‌ها در نقش `microbiome_counts.csv` قرار می‌گیرد و حجم‌های مغزی در نقش `brain_features.csv`. بخش‌های sCCA، تحلیل میانجی و رگرسیون روی نمرهٔ شناختی همه قابل اجرا هستند.

### D2 — Khula: میکروبیوم + **EEG** در ۱۹۴ نوزاد (مناسب تخصص EEG تو)

| | |
|---|---|
| **مقاله** | Bonham, Margolis et al., 2025, *mBio*، با عنوان *Co-development of gut microbial metabolism and visual neural circuitry over human infancy* |
| **افراد** | ۱۹۴ نوزاد در کیپ‌تاون، با ۳ ویزیت در حدود ۴، ۹ و ۱۴ ماهگی |
| **داده** | ژن‌های میکروبی از متاژنوم شاتگان (`bonham_margolis_mbio_allunirefs.csv.gz`)، و **پتانسیل برانگیختهٔ بینایی (VEP)** شامل دامنه و تأخیر مؤلفه‌های N1، P1 و N2 در فایل `allmeta.csv` |
| **دسترسی** | دانلود آزاد از Dryad، حدود ۵۱۱ مگابایت: [doi:10.5061/dryad.rn8pk0pq2](https://datadryad.org/dataset/doi:10.5061/dryad.rn8pk0pq2). کد: [Zenodo 10.5281/zenodo.15723691](https://doi.org/10.5281/zenodo.15723691). مقاله: [PMC12345167](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12345167/) |

**ارتباط با اوتیسم:** این طرح فرضیه‌محور است. ماژول‌های ژنی روده–مغز که در اوتیسم مطرح شده‌اند را انتخاب کن، مثل سنتز و تخریب GABA، مسیر تریپتوفان و کینورنین، تولید propionate و p-cresol، و متابولیسم گوگرد. بعد بررسی کن آیا فراوانی این ماژول‌ها در ۴ ماهگی، تکامل VEP در ۹ تا ۱۴ ماهگی را پیش‌بینی می‌کند. پردازش بینایی اولیه یکی از نشانگرهای زودهنگام مطرح در اوتیسم است.
محدودیت: فقط ویژگی‌های مشتق‌شدهٔ VEP منتشر شده و EEG خام در دسترس نیست. پس source localization ممکن نیست. پروفایل گونه‌ها هم در این فایل‌ها نیست و فقط ژن‌ها هستند.

### D3 — Mendelian randomization: روده → مغز → اوتیسم با آمار خلاصهٔ GWAS

| منبع | محتوا | دسترسی |
|---|---|---|
| **MiBioGen** (Kurilshikov et al., 2021, *Nat Genet*) | GWAS فراوانی taxaهای روده در ۱۸٬۳۴۰ نفر | دانلود از سایت MiBioGen. [مقالهٔ کنسرسیوم](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5992867/) |
| **BIG40** (Smith et al., 2021, *Nat Neurosci*) | GWAS روی **۳۹۳۵ فنوتیپ تصویربرداری مغز** در UK Biobank (کوهورت کشف حدود ۲۲ هزار نفر) | دانلود آزاد: [open.win.ox.ac.uk/ukbiobank/big40](https://open.win.ox.ac.uk/ukbiobank/big40) و EBI GWAS Catalog با شناسه‌های GCST90002426 تا GCST90006360 |
| **iPSYCH-PGC ASD** (Grove et al., 2019, *Nat Genet*) | GWAS تشخیص اوتیسم با ۱۸٬۳۸۱ بیمار و ۲۷٬۹۶۹ کنترل | دانلود از [PGC](https://pgc.unc.edu/for-researchers/download-results/) بعد از پذیرفتن شرایط استفاده. [شرایط iPSYCH](https://ipsych.au.dk/en/research/downloads) |

**ارتباط با اوتیسم:** این تنها طرح کاملاً باز است که پیامد آن **خودِ تشخیص اوتیسم** است. طرح پیشنهادی MR دومرحله‌ای (mediation MR) است: taxon → فنوتیپ مغزی (IDP) → ASD. برای اطمینان از نتایج، تحلیل حساسیت (MR-Egger، weighted median، MR-PRESSO) و colocalization هم لازم است.
در این طرح داده‌ی فردی وجود ندارد و EEG هم نیست. مقاله‌های MR مشابه منتشر شده‌اند، از جمله [NeuroImage 2025](https://www.sciencedirect.com/science/article/pii/S1053811925002952) و [MR روده–اوتیسم ۲۰۲۳](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10753022/). پس نوآوری باید از جای دیگری بیاید، مثلاً از mediation چندمتغیره روی IDPهای مشخص یا استفاده از GWASهای جدیدتر میکروبیوم.

### D4 — مدل موشی اوتیسم: روده + مغز + رفتار در **همان حیوانات**

| | |
|---|---|
| **مقاله** | Sharon et al., 2019, *Cell*، با عنوان *Human gut microbiota from ASD promote behavioral symptoms in mice* |
| **طرح** | میکروبیوتای کودکان ASD و TD به موش‌های بدون میکروب پیوند زده شد |
| **داده** | دادهٔ سطح شکل‌ها، شامل رفتار و 16S: [Mendeley Data](https://data.mendeley.com/datasets/ngzmj4zkms). متابولومیکس پلاسما و کولون: MetaboLights [MTBLS726](https://www.elucidata.io/datasets/dataset-of-the-week-mtbls726). **RNA-seq مغز**: GEO [GSE109827](https://pluto.bio/product/curated-data-sets/explore/experiment/gse109827-gut-microbiota-from-human-autism-spectrum-disorder-induces-behavioral-deficits-in-mice) |

**ارتباط با اوتیسم:** این طرح مستقیماً به اوتیسم مربوط است و روده، متابولیت، مغز و رفتار را در همان حیوانات دارد. مغز در این‌جا ترانسکریپتوم است و تصویربرداری یا EEG ندارد.

### داده‌های باز اوتیسم برای ساختن «امضای میکروبی» (ورودی طرح D1)

| دیتاست | محتوا |
|---|---|
| **Kang 2017 MTT** ([figshare](https://springernature.figshare.com/collections/Microbiota_Transfer_Therapy_alters_gut_ecosystem_and_improves_gastrointestinal_and_autism_symptoms_an_open-label_study/3674866)) | 16S، متابولیت‌ها، CARS، SRS، ABC و GSRS |
| **GEO GSE113540** ([OmicsDI](https://omicsdi.org/dataset/geo/GSE113540)) | متاژنوم شاتگان مدفوع از ۳۰ کودک ASD و ۳۰ کنترل. فرمت آن با RESONANCE سازگارتر است |
| **PRJNA952791**، **GSE113690** و کوهورت‌های فراتحلیل Sci Rep 2022 | 16S در چند کوهورت؛ برای اعتبارسنجی بین‌کوهورتی امضا (leave-one-cohort-out) |

> نکتهٔ هماهنگ‌سازی: داده‌های 16S را فقط در **سطح جنس** با داده‌های شاتگان مقایسه کن. داده‌های شاتگان را با یک ابزار واحد پردازش کن، مثلاً MetaPhlAn 4 از روی ریدهای خام SRA. امضا را هم فقط روی داده‌ای بساز که سن کودکانش با کودکان هدف هم‌خوان باشد، چون میکروبیوم نوزاد با میکروبیوم کودک ۶ ساله خیلی فرق دارد.

---

## سطح ۱ — مغز + روده در همان کودکان (بهترین تطابق با هدف تو)

| مطالعه | افراد | داده‌ها | دسترسی |
|---|---|---|---|
| **Mi et al., 2026, *Cell Reports Medicine*** — *An integrative multi-omics approach identifies microbiome alterations linked to pathological and behavioral features in ASD* | ۳۲۶ کودک ASD + ۱۶۹ کودک TD، سن ۰ تا ۱۰ سال | **MRI ساختاری** + متاژنوم روده + متابولوم پلاسما + شدت علائم | بخش Data availability می‌گوید ریدهای متاژنومی در یک پایگاه ملی چین ثبت شده‌اند (شمارهٔ دقیق را از مقاله بردار)؛ دادهٔ تصویربرداری احتمالاً **با درخواست از lead contact**. [مقاله](https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(26)00072-8) · [PubMed 41806837](https://pubmed.ncbi.nlm.nih.gov/41806837/) · DOI: 10.1016/j.xcrm.2026.102655 |
| **Aziz-Zadeh et al., 2025, *Nature Communications*** — *Relationships between brain activity, tryptophan-related gut metabolites, and autism symptomatology* | ۴۳ ASD + ۴۱ نوروتیپیک، ۸ تا ۱۷ سال | **fMRI تکلیف‌محور** (تکالیف اجتماعی–هیجانی و حسی) + متابولومیکس مدفوع + ارزیابی رفتاری | تأمین مالی NIH (NICHD R01HD079432) ⇒ ممکن است در **NIMH Data Archive (NDA)** هم ثبت شده باشد؛ بخش Data availability مقاله را چک کن. [مقاله](https://www.nature.com/articles/s41467-025-58459-1) · [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11302680/) · [صفحهٔ پروژه در USC](https://chan.usc.edu/research/projects/the-relationship-between-brain-functioning-behavior-and-microbiota-in-autism-spectrum-disorder) |
| بیمارستان دهم شانگهای — FMT و rs-fMRI (ALFF و ReHo) | ۸۰ کودک ASD، ۶ ماه پس از پیوند میکروبیوتای مدفوع | rs-fMRI قبل/بعد از FMT | فقط چکیدهٔ کنفرانس؛ برای دادهٔ کامل باید با نویسندگان تماس گرفت. [ISMRM 2023](https://archive.ismrm.org/2023/5303.html) · [ISMRM 2024](https://archive.ismrm.org/2024/4403.html) |

**یافته‌های کلیدی برای فرضیه‌سازی:** در Mi و همکاران، فراوانی *Clostridioides difficile*
قوی‌ترین پیش‌بین شدت علائم و تغییرات ساختاری مغز بود و تفاوت‌ها با افزایش سن کم می‌شد.
در Aziz-Zadeh و همکاران، متابولیت‌های مسیر تریپتوفان (مثل kynurenate، کمتر در ASD) با
فعالیت **اینسولا** و **سینگولیت** و شدت علائم مرتبط بودند ⇒ مسیر
«میکروب → متابولیت تریپتوفان → شبکهٔ salience → علائم» فرضیهٔ خوبی برای مدل میانجی است.

> متن نمونهٔ ایمیل درخواست داده در انتهای همین فایل آمده است.

---

## سطح ۲ — میکروبیوم + آزمایش‌ها/مقیاس‌های بالینی (داده‌های باز)

| دیتاست | افراد | داده‌ها | دسترسی |
|---|---|---|---|
| **Kang et al., 2017, *Microbiome* — Microbiota Transfer Therapy (MTT)** | ۱۸ ASD (+۲۰ TD در مطالعهٔ متابولومیکس پیگیری) | 16S، بعداً متاژنوم شاتگان و متابولومیکس مدفوع (۶۶۹ ترکیب)، **CARS، SRS، ABC، PGI-II، GSRS**؛ طولی (قبل، حین، بعد از درمان و پیگیری ۲ ساله) | [figshare](https://springernature.figshare.com/collections/Microbiota_Transfer_Therapy_alters_gut_ecosystem_and_improves_gastrointestinal_and_autism_symptoms_an_open-label_study/3674866) (DOI 10.6084/m9.figshare.c.3674866.v1) · کد: `github.com/gregcaporaso/autism-fmtl` · [متاژنوم](https://pmc.ncbi.nlm.nih.gov/articles/PMC9654974/) · [ادغام چنداُمیکسی MTT، bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681677.full.pdf) |
| **Morton et al., 2023, *Nature Neuroscience*** — *Multi-level analysis of the gut–brain axis…* | ۱۰ دیتاست میکروبیوم + ۱۵ دیتاست دیگر | رژیم غذایی، متابولومیکس، سایتوکاین‌ها، **بیان ژن مغز انسان**؛ الگوی مشترک در کوهورت‌های هم‌سن/هم‌جنس دیده شد ولی در کوهورت‌های خواهر–برادری نه | منابع همه عمومی‌اند (بخش Data/Code availability مقاله). DOI: 10.1038/s41593-023-01361-0 · [متن](https://escholarship.org/content/qt06b011nh/qt06b011nh.pdf) |
| **Su et al., 2024, *Nature Microbiology*** — نشانگرهای چندقلمرویی (باکتری، قارچ، ویروس، آرکیا) | ۱۶۲۷ کودک | متاژنوم شاتگان + اعتبارسنجی خارجی روی ۶ مطالعهٔ منتشرشده | [مقاله](https://www.nature.com/articles/s41564-024-01739-1) — شمارهٔ دسترسی در Data availability |
| **Dan et al., 2020, *Gut Microbes*** | ۱۴۳ کودک ۲ تا ۱۳ ساله | 16S + متابولیت‌های مدفوع | [DOAJ](https://doaj.org/article/bde4e2e2c0814bd895437c65e8898ef5) |
| **Zhong et al., 2026, *Frontiers in Microbiology*** | کودکان ASD و TD (تیانجین) | **توالی‌یابی کل اگزوم** + 16S + متابولومیکس پلاسما؛ واریانت‌های ژن‌های MUC | [مقاله](https://www.frontiersin.org/articles/10.3389/fmicb.2026.1766850/full) · [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12913364/) |
| **BioProject PRJNA952791** (دانشگاه یوتا) | ASD + اعضای نوروتیپیک خانواده | 16S مدفوع + موش‌های gnotobiotic کلونیزه‌شده با همین میکروبیوتا (۱۶۷ آزمایش SRA) | [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/952791) |
| **GEO GSE113690** | ASD در برابر سالم | 16S rRNA | [OmicsDI](https://omicsdi.org/dataset/geo/GSE113690) |
| **فراتحلیل Sci Rep 2022** | ۱۳ دیتاست از ۱۰ کوهورت | 16S و متاژنوم عمومی (SRA) — بهترین نقطهٔ شروع برای یافتن accessionها و اعتبارسنجی بین‌کوهورتی | [مقاله](https://www.nature.com/articles/s41598-022-21327-9) |
| **یادگیری ماشین روی 16S عمومی (۲۰۲۵)** | ۶۹۲ نمونه (چین، آمریکا، اکوادور) | 16S تجمیعی | [PMC12765537](https://pmc.ncbi.nlm.nih.gov/articles/PMC12765537/) |
| **Yap et al., 2021, *Cell*** — Australian Autism Biobank + QTAB | ۲۴۷ کودک | متاژنوم + رژیم غذایی + فنوتیپ؛ نشان داد **رژیم غذایی** ارتباط اوتیسم–میکروبیوم را میانجی‌گری می‌کند | با درخواست از Australian Autism Biobank. [مقاله](https://www.sciencedirect.com/science/article/pii/S0092867421012319) |

---

## سطح ۳ — فقط مغز (برای ساخت و اعتبارسنجی شاخهٔ مغزی مدل)

| دیتاست | داده‌ها | دسترسی |
|---|---|---|
| **Healthy Brain Network (HBN)**، Child Mind Institute | کودکان و نوجوانان ۵ تا ۲۱ سال (حدود ۵۰۰۰ نفر)؛ MRI ساختاری و عملکردی، **EEG پرچگال ۱۲۸ کاناله** + eye-tracking، ارزیابی‌های روان‌پزشکی؛ تعداد قابل‌توجهی تشخیص اوتیسم دارند (فیلدهای تشخیص را در فایل‌های فنوتیپ هر release چک کن؛ دادهٔ فنوتیپ معمولاً به توافق‌نامهٔ استفاده نیاز دارد) | [EEG](https://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/downloads/downloads_EEG_R2_1.html) · [MRI](https://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/downloads/downloads_MRI_R3.html) · AWS S3 بی‌نام: `s3://fcp-indi/data/Projects/HBN` · [HBN-EEG در قالب BIDS](https://www.biorxiv.org/content/10.1101/2024.10.03.615261v2) · [EEGDash](https://eegdash-readme.static.hf.space/index.html) |
| **ABIDE I / ABIDE II** | rs-fMRI + T1 از ده‌ها سایت، ASD و کنترل؛ نسخهٔ پیش‌پردازش‌شدهٔ PCP مستقیماً با nilearn دانلود می‌شود | `gutbrain.fmri.fetch_abide_features()` در همین مخزن |

**کاربرد:** با HBN-EEG شاخهٔ EEG مدل را (همان ویژگی‌های `gutbrain/eeg.py`) روی صدها کودک
آموزش و اعتبارسنجی کن؛ بعد همان پایپ‌لاین را روی دادهٔ زوج خودت (مسیر C) یا دادهٔ درخواستی
(مسیر A) اجرا کن. این کار از overfitting روی نمونهٔ کوچک دادهٔ زوج جلوگیری می‌کند.

---

## سطح ۴ — کوهورت‌های آینده‌نگر نوزادان (برای کارهای بعدی)

| کوهورت | توضیح |
|---|---|
| **GEMMA** (Genome, Environment, Microbiome and Metabolome in Autism) | نوزادان ۰ تا ۶ ماهه با خواهر/برادر مبتلا؛ خون، مدفوع، ادرار، بزاق تا ۳ سالگی. [MGH](https://www.massgeneral.org/children/research/genome-environment-microbiome-metabolome-in-autism-study) · [NCT04271774](https://clinicaltrials.gov/study/NCT04271774) |
| **EASE** (سوئد) | نوزادان با احتمال بالای ASD از ۵ تا ۳۶ ماهگی. [Transl Psychiatry 2023](https://www.nature.com/articles/s41398-023-02556-6) |
| **MARBLES** | شواهد مختلط برای ارتباط میکروبیوم نوزادی و ASD بعدی. [Autism Research 2026](https://onlinelibrary.wiley.com/doi/10.1002/aur.70207) |
| **ECHO** | میکروبیوم ۱ سالگی و رفتار اجتماعی (SRS-2) در ۳ سالگی. [PMC11101342](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11101342/) |
| **NDA collection 2557** | مطالعهٔ آینده‌نگر خواهر/برادرهای پرخطر و میکروبیوم روده. [NDA](https://nda.nih.gov/edit_collection.html?id=2557) |

---

## سه مسیر پیشنهادی

### مسیر A — درخواست دادهٔ زوج (۲ تا ۶ ماه زمان)
به lead contact مقالهٔ Mi و همکاران (۲۰۲۶) و نویسندگان Aziz-Zadeh و همکاران (۲۰۲۵) ایمیل بزن.
معمولاً یک **Data Use Agreement** و تأیید کمیتهٔ اخلاق دانشگاهت لازم است. برای NDA باید
از طریق یک مؤسسهٔ دارای حساب NIH درخواست داد.

### مسیر B — شروع فوری با داده‌های باز
- شاخهٔ روده: Kang 2017 (تنها دیتاست باز با **میکروبیوم + متابولیت + CARS/SRS/ABC/GSRS** در همان کودکان)
  + کوهورت‌های 16S عمومی برای اعتبارسنجی **leave-one-cohort-out**.
- شاخهٔ مغز: HBN-EEG و ABIDE.
- پل بین دو شاخه: فقط در سطح **فرضیه** (مثلاً مسیر تریپتوفان/کینورنین، SCFAها، p-cresol)،
  نه در سطح فرد.

### مسیر C — جمع‌آوری دادهٔ خودت (پیشنهاد اصلی برای پایان‌نامه/مقاله)
پروتکل پیشنهادی در [`model.md`](model.md#پروتکل-جمعآوری-داده) آمده است.

---

## نمونهٔ ایمیل درخواست داده

```
Subject: Data access request – gut microbiome and neuroimaging data in children with ASD

Dear Dr. <Name>,

I am a biomedical engineering researcher at <University>, working on EEG/neuroimaging
analysis of the gut–brain axis in autism. I read your article "<title>" (<journal>, <year>)
with great interest.

I would like to request access to the de-identified data from this study (gut microbiome
profiles, <neuroimaging/fMRI-derived> features, metabolomics and clinical scores) to test
a multimodal integration framework (sparse CCA, mediation and nested cross-validated
fusion models). The analysis code is open source: <link to this repository>.

I am happy to sign a data use agreement, obtain approval from our ethics committee,
cite your work, and share all results with your team before submission. Collaboration
or co-authorship can be discussed as appropriate.

Thank you for considering this request.

Best regards,
<Name>, <Position>, <Affiliation>, <Email>
```
