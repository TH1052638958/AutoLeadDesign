import openai
import json



# model_engine = 'gpt-3.5-turbo-instruct' # You can choose a different model if desired
#
# prompt = f'Generate a novel valid molecule  and do not generate any English text'
# response = openai.Completion.create(
#             engine=model_engine,
#             prompt=prompt,
#             max_tokens=60,
#             temperature=0.7,
#             n=1,
#             stop=None,
#             timeout=20
#         )
# new_mol = response.choices[0].text.strip()
#
# print("new molecules", new_mol)

import os
from openai import OpenAI
from openai.types.chat import completion_create_params
os.environ["OPENAI_BASE_URL"]="https://api.deepseek.com"
os.environ["OPENAI_API_KEY"]="xxx"
model_engine = 'deepseek-v4-flash' # You can choose a different model if desired
prompt = f'Generate a novel valid drug-like molecule SMILES which contains one fragment of [ N1CC[C@@]2(C[C@H](N)C3=C(C=CC=C3)O2)C1 , C1NC[C@@H]2CC[C@]3(N=C(CC)NC3=O)[C@@H]12,N1CC2=CC(C)=CC=C2C1C1=CC=CC=C1 ] at least and do not generate any English text.'
client = OpenAI()
completion = client.chat.completions.create(
                    model=model_engine,
                    messages=[
                        {"role": "user", "content": prompt}],
                    #n=1,
                    #max_tokens=60,
                    temperature=1.5,
                    extra_body={"thinking": {"type": "disabled"}},
                    #stop="!",
                    #user="user"
                    stream=False
                )
new_mol = json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
print(new_mol)



