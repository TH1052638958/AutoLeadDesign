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
os.environ["OPENAI_BASE_URL"]="https://www.blueshirtmap.com/v1"
os.environ["OPENAI_API_KEY"]="xxx"
model_engine = 'gpt-4o' # You can choose a different model if desired
prompt = f'Generate a novel valid molecule SMILES and do not generate any English text'

client = OpenAI()

completion = client.chat.completions.create(
                model=model_engine,
                messages=[
                    {"role": "user", "content": prompt}],
                n=1,
                max_tokens=60,
                temperature=0.5,
                stop="!",
                user="user"
            )
print(completion)
new_mol=json.loads(completion.model_dump_json())["choices"][0]["message"]["content"]
print(new_mol)


