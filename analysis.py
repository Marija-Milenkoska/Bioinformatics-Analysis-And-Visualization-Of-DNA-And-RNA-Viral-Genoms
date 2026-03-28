
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


data = pd.read_csv("virus_dataset_project.csv")

#analysis part

print("HEAD:")
print(data.head(), "\n")

print("INFO:")
print(data.info(), "\n")

print("MISSING VALUES:")
print(data.isna().sum(), "\n")

print("DUPLICATES:")
print(data.duplicated().sum(), "\n")

data["Treatment"] = data["Treatment"].fillna("None")
data["Vaccine"] = data["Vaccine"].fillna("None")

print("GENOME TYPE COUNTS:")
print(data["Genome_type"].value_counts(), "\n")

print("HOST FREQUENCY:")
print(data["Host"].value_counts(), "\n")

print("STRAND TYPE FREQUENCY:")
print(data["Strand_type"].value_counts(), "\n")

print("DESCRIPTIVE STATISTICS:")
print(data[["Genome_length", "GC_content"]].describe(), "\n")

print("AVERAGE GENOME LENGTH BY GENOME TYPE:")
print(data.groupby("Genome_type")["Genome_length"].mean(), "\n")

print("AVERAGE GC CONTENT BY GENOME TYPE:")
print(data.groupby("Genome_type")["GC_content"].mean(), "\n")

print("HOST DISTRIBUTION BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Host"]), "\n")

print("VACCINE AVAILABILITY BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Vaccine"]), "\n")

print("TREATMENT AVAILABILITY BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Treatment"]), "\n")

print("LONGEST GENOMES:")
print(data.nlargest(5, "Genome_length")[["Virus_name", "Genome_type", "Genome_length"]], "\n")

print("SHORTEST GENOMES:")
print(data.nsmallest(5, "Genome_length")[["Virus_name", "Genome_type", "Genome_length"]], "\n")

print("FAMILIES PER GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Family"]), "\n")

print("DISEASE DISTRIBUTION:")
print(data["Disease"].value_counts(), "\n")

print("DISEASE DISTRIBUTION BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Disease"]), "\n")

print("TRANSMISSION DISTRIBUTION:")
print(data["Transmission"].value_counts(), "\n")

print("TRANSMISSION BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Transmission"]), "\n")

print("TARGET ORGAN DISTRIBUTION:")
print(data["Target_organ"].value_counts(), "\n")

print("TARGET ORGAN BY GENOME TYPE:")
print(pd.crosstab(data["Genome_type"], data["Target_organ"]), "\n")

#visualization part

sns.set(style="whitegrid")

plt.figure(figsize=(8, 5))
sns.countplot(data=data, x="Genome_type", hue="Genome_type", legend=False)
plt.title("Број на DNA и RNA вируси")
plt.xlabel("Тип на геном")
plt.ylabel("Број на вируси")
plt.tight_layout()
plt.savefig("genome_type_count.png")
plt.close()


plt.figure(figsize=(12, 6))
sns.barplot(
    data=data.sort_values("Genome_length", ascending=False),
    x="Virus_name",
    y="Genome_length",
    hue="Genome_type"
)
plt.title("Должина на геном кај различни вируси")
plt.xlabel("Вирус")
plt.ylabel("Должина на геном")
plt.xticks(rotation=75, ha="right")
plt.tight_layout()
plt.savefig("genome_length_by_virus.png")
plt.close()

avg_length = data.groupby("Genome_type", as_index=False)["Genome_length"].mean()
plt.figure(figsize=(8, 5))
sns.barplot(data=avg_length, x="Genome_type", y="Genome_length", hue="Genome_type", legend=False)
plt.title("Просечна должина на геном кај DNA и RNA вируси")
plt.xlabel("Тип на геном")
plt.ylabel("Просечна должина на геном")
plt.tight_layout()
plt.savefig("average_genome_length.png")
plt.close()

plt.figure(figsize=(12, 6))
sns.barplot(
    data=data.sort_values("GC_content", ascending=False),
    x="Virus_name",
    y="GC_content",
    hue="Genome_type"
)
plt.title("GC содржина кај различни вируси")
plt.xlabel("Вирус")
plt.ylabel("GC содржина (%)")
plt.xticks(rotation=75, ha="right")
plt.tight_layout()
plt.savefig("gc_content_by_virus.png")
plt.close()

avg_gc = data.groupby("Genome_type", as_index=False)["GC_content"].mean()
plt.figure(figsize=(8, 5))
sns.barplot(data=avg_gc, x="Genome_type", y="GC_content", hue="Genome_type", legend=False)
plt.title("Просечна GC содржина кај DNA и RNA вируси")
plt.xlabel("Тип на геном")
plt.ylabel("Просечна GC содржина (%)")
plt.tight_layout()
plt.savefig("average_gc_content.png")
plt.close()

host_counts = pd.crosstab(data["Genome_type"], data["Host"])
host_counts.plot(kind="bar", figsize=(10, 6))
plt.title("Распределба на домаќини според тип на геном")
plt.xlabel("Тип на геном")
plt.ylabel("Број на вируси")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("host_distribution.png")
plt.close()

plt.figure(figsize=(9, 5))
sns.countplot(data=data, x="Strand_type", hue="Genome_type")
plt.title("Распределба на типови на нишка")
plt.xlabel("Тип на нишка")
plt.ylabel("Број на вируси")
plt.tight_layout()
plt.savefig("strand_type_distribution.png")
plt.close()

vaccine_counts = data["Vaccine"].value_counts()
plt.figure(figsize=(7, 7))
plt.pie(vaccine_counts, labels=vaccine_counts.index, autopct="%1.1f%%", startangle=90)
plt.title("Достапност на вакцина")
plt.tight_layout()
plt.savefig("vaccine_availability.png")
plt.close()

plt.figure(figsize=(8, 5))
sns.boxplot(data=data, x="Genome_type", y="Genome_length", hue="Genome_type", legend=False)
plt.title("Распределба на должина на геном според тип")
plt.xlabel("Тип на геном")
plt.ylabel("Должина на геном")
plt.tight_layout()
plt.savefig("genome_length_boxplot.png")
plt.close()

plt.figure(figsize=(8, 5))
sns.boxplot(data=data, x="Genome_type", y="GC_content", hue="Genome_type", legend=False)
plt.title("Распределба на GC содржина според тип")
plt.xlabel("Тип на геном")
plt.ylabel("GC содржина (%)")
plt.tight_layout()
plt.savefig("gc_content_boxplot.png")
plt.close()

plt.figure(figsize=(10, 5))
data["Disease"].value_counts().plot(kind="bar")
plt.title("Распределба на вируси според болест")
plt.xlabel("Болест")
plt.ylabel("Број на вируси")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("disease_distribution.png")
plt.close()

pd.crosstab(data["Genome_type"], data["Disease"]).plot(kind="bar", figsize=(12,6))
plt.title("Распределба на болести според тип на геном")
plt.xlabel("Тип на геном")
plt.ylabel("Број на вируси")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("disease_by_genome.png")
plt.close()

plt.figure(figsize=(10, 5))
data["Transmission"].value_counts().plot(kind="bar")
plt.title("Распределба на вируси според начин на пренос")
plt.xlabel("Начин на пренос")
plt.ylabel("Број на вируси")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("transmission_distribution.png")
plt.close()

pd.crosstab(data["Genome_type"], data["Transmission"]).plot(kind="bar", figsize=(12,6))
plt.title("Распределба на начин на пренос според тип на геном")
plt.xlabel("Тип на геном")
plt.ylabel("Број на вируси")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("transmission_by_genome.png")
plt.close()

plt.figure(figsize=(10, 5))
data["Target_organ"].value_counts().plot(kind="bar")
plt.title("Распределба на вируси според целен орган")
plt.xlabel("Целен орган")
plt.ylabel("Број на вируси")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("target_organ_distribution.png")
plt.close()

pd.crosstab(data["Genome_type"], data["Target_organ"]).plot(kind="bar", figsize=(12,6))
plt.title("Распределба на целен орган според тип на геном")
plt.xlabel("Тип на геном")
plt.ylabel("Број на вируси")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("target_organ_by_genome.png")
plt.close()

plt.figure(figsize=(12, 6))
sorted_data = data.sort_values("Genome_length")
plt.plot(sorted_data["Virus_name"], sorted_data["Genome_length"], marker='o')
plt.title("Линеарен приказ на должината на геномот кај вируси")
plt.xlabel("Вирус")
plt.ylabel("Должина на геном")
plt.xticks(rotation=75, ha="right")
plt.tight_layout()
plt.savefig("genome_length_line.png")
plt.close()

plt.figure(figsize=(12, 6))
sorted_gc = data.sort_values("GC_content")
plt.plot(sorted_gc["Virus_name"], sorted_gc["GC_content"], marker='o')
plt.title("Линеарен приказ на GC содржината кај вируси")
plt.xlabel("Вирус")
plt.ylabel("GC содржина (%)")
plt.xticks(rotation=75, ha="right")
plt.tight_layout()
plt.savefig("gc_content_line.png")
plt.close()