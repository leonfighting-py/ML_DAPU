import csv
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / 'DA-PU数据整合-new.csv'
OUTPUT_FILE = BASE_DIR / 'DA-PU数据整合-new_furan_cleaned.csv'


def normalize_text(value):
    return (value or '').strip()


def classify_furan_diol(component_name):
    name = normalize_text(component_name)
    name_lower = name.lower()

    if any(keyword in name for keyword in ['呋喃二甲醇', '2,5呋喃二甲醇']) or 'fdm' in name_lower:
        return '呋喃二甲醇'

    if any(keyword in name for keyword in ['糠醇']) or 'kc' in name_lower:
        return '糠醇'

    if any(keyword in name for keyword in ['糠胺']) or any(keyword in name_lower for keyword in ['ka', 'fge', 'fdo', 'fpdp']):
        return '糠胺'

    return ''


def main():
    with INPUT_FILE.open('r', encoding='utf-8-sig', newline='') as source:
        reader = csv.reader(source)
        rows = list(reader)

    if not rows:
        raise ValueError('输入 CSV 为空，无法执行清洗。')

    header = [str(cell).strip() for cell in rows[0]]
    try:
        component_idx = header.index('组分名称component_name')
        group_idx = header.index('组分角色group_type（Isocyanate Hydroxyl Amine')
    except ValueError as exc:
        raise ValueError('找不到所需列：`组分名称component_name` 或 `组分角色group_type（Isocyanate Hydroxyl Amine`。') from exc

    output_header = header + ['furan_diol_clean']
    cleaned_rows = [output_header]
    summary = {'呋喃二甲醇': 0, '糠醇': 0, '糠胺': 0, '': 0}

    for row in rows[1:]:
        padded_row = list(row)
        if len(padded_row) < len(header):
            padded_row.extend([''] * (len(header) - len(padded_row)))

        group_type = normalize_text(padded_row[group_idx]).lower()
        component_name = padded_row[component_idx]

        cleaned_label = ''
        if 'furan_diol' in group_type:
            cleaned_label = classify_furan_diol(component_name)

        summary[cleaned_label] = summary.get(cleaned_label, 0) + 1
        cleaned_rows.append(padded_row + [cleaned_label])

    with OUTPUT_FILE.open('w', encoding='utf-8-sig', newline='') as target:
        writer = csv.writer(target)
        writer.writerows(cleaned_rows)

    print(f'已生成清洗文件: {OUTPUT_FILE.name}')
    print('furan_diol_clean 分布:')
    for label in ['呋喃二甲醇', '糠醇', '糠胺', '']:
        display = label if label else '(空)'
        print(f'  {display}: {summary.get(label, 0)}')


if __name__ == '__main__':
    main()